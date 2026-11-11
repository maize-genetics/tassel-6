//! Times the shared WGSL Gram and a platform BLAS on the four published shapes.
//!
//! `cargo run --release` writes `target/gpu-gram-report.json`. `--blas-only` is the one-thread
//! child, so Accelerate or OpenBLAS sees its thread-count variable before it initializes.

mod blas;
mod gpu;
mod report;

use std::path::{Path, PathBuf};
use std::process::Command;
use std::time::Instant;

use gpu_gram::{
    dosage_checksum, dosage_gram, generate_dosages, PARK_MILLER_SEED, PANEL_COLUMNS, PREVIEW, SHAPES,
};

use crate::blas::BlasInput;
use crate::report::{AdapterReport, BlasRow, GpuRow, HostReport, Report, SizeReport};

fn main() {
    let mut args = std::env::args().skip(1);
    if args.next().as_deref() == Some("--blas-only") {
        let path = args.next().expect("--blas-only needs an output path");
        write_one_thread(Path::new(&path));
        return;
    }

    let report_path = Path::new(env!("CARGO_MANIFEST_DIR")).join("target/gpu-gram-report.json");
    std::fs::create_dir_all(report_path.parent().unwrap()).expect("target directory");
    let mut report = Report {
        host: host_report(),
        adapter: None,
        gpu_error: None,
        setup_seconds: None,
        blas_library: None,
        sizes: Vec::new(),
    };
    let gpu = match gpu::open() {
        Ok(gpu) => {
            if gpu.fallback {
                eprintln!(
                    "wgpu adapter {} ({}) is a software or virtual device; timings are labelled fallback",
                    gpu.adapter_name, gpu.device_type
                );
            } else {
                eprintln!(
                    "wgpu adapter {} backend {} type {}",
                    gpu.adapter_name, gpu.backend, gpu.device_type
                );
            }
            report.setup_seconds = Some(gpu.setup_seconds);
            report.adapter = Some(AdapterReport {
                name: gpu.adapter_name.clone(),
                vendor: gpu.vendor,
                device: gpu.device_id,
                device_type: gpu.device_type.clone(),
                backend: gpu.backend.clone(),
                driver: gpu.driver.clone(),
                driver_info: gpu.driver_info.clone(),
                fallback: gpu.fallback,
                limits: gpu.limits.clone(),
            });
            Some(gpu)
        }
        Err(error) => {
            eprintln!("wgpu unavailable: {error}");
            report.gpu_error = Some(error);
            None
        }
    };

    let mut prepared = Vec::new();
    for shape in SHAPES {
        eprint!("{}x{}: generating dosages... ", shape.taxa, shape.sites);
        let codes = generate_dosages(shape.taxa, shape.sites, PARK_MILLER_SEED);
        let checksum = dosage_checksum(&codes);
        assert_eq!(checksum, shape.checksum, "{}x{} checksum", shape.taxa, shape.sites);
        assert_eq!(&codes[..PREVIEW.len()], &PREVIEW);
        eprintln!("checksum {checksum}");
        eprint!("{}x{}: f64 reference... ", shape.taxa, shape.sites);
        let started = Instant::now();
        let reference = dosage_gram(&codes, shape.taxa, shape.sites, PANEL_COLUMNS);
        let reference_seconds = started.elapsed().as_secs_f64();
        eprintln!("{reference_seconds:.3} s");
        cache_reference(shape.taxa, shape.sites, &reference);
        let item = Prepared {
            taxa: shape.taxa,
            sites: shape.sites,
            checksum,
            codes,
            reference,
            reference_seconds,
        };
        let rows = if let Some(gpu) = &gpu {
            eprint!("{}x{}: wgpu... ", item.taxa, item.sites);
            match gpu::time_shape(gpu, &item.codes, item.taxa, item.sites, &item.reference) {
                Ok(timings) => {
                    for timing in &timings {
                        eprintln!(
                            "{} median {:.6} s max_abs {:.3e} bit_exact {}",
                            timing.shader,
                            timing.phases.end_to_end.median_seconds,
                            timing.error.max_abs,
                            timing.error.bit_exact,
                        );
                    }
                    gpu_rows(&timings)
                }
                Err(error) => {
                    eprintln!("failed: {error}");
                    report.gpu_error = Some(error);
                    Vec::new()
                }
            }
        } else {
            Vec::new()
        };
        report.sizes.push(size_report(&item, rows, Vec::new()));
        prepared.push(item);
    }

    let blas_name = blas::backend_label();
    let thread_env = blas::thread_env_var();
    if let Ok(previous) = std::env::var(thread_env) {
        eprintln!("unsetting {thread_env}={previous} for the default {blas_name} run");
        std::env::remove_var(thread_env);
    }
    let library_line = blas::library_line();
    eprintln!("{library_line}");
    report.blas_library = Some(library_line);
    eprintln!("{blas_name}, default thread pool");
    for (index, item) in prepared.iter().enumerate() {
        let outcomes = blas::time_shape(
            &BlasInput {
                taxa: item.taxa,
                sites: item.sites,
                codes: &item.codes,
                reference: &item.reference,
            },
            "default",
        );
        for outcome in &outcomes {
            eprintln!(
                "{}x{} {} median {:.6} s max_abs {:.3e}",
                item.taxa, item.sites, outcome.kernel, outcome.syrk.median_seconds, outcome.max_abs
            );
        }
        report.sizes[index].blas.extend(outcomes.into_iter().map(blas_row));
    }
    drop(prepared);

    let one_thread_path = report_path.with_file_name("blas-1thread.txt");
    eprintln!("{blas_name}, {thread_env}=1");
    let status = Command::new(std::env::current_exe().expect("current exe"))
        .arg("--blas-only")
        .arg(&one_thread_path)
        .env(thread_env, "1")
        .status()
        .expect("spawn one-thread blas");
    assert!(status.success(), "one-thread blas exited {status}");
    let text = std::fs::read_to_string(&one_thread_path).expect("one-thread blas output");
    for (taxa, sites, outcome) in blas::parse_lines(&text) {
        eprintln!(
            "{taxa}x{sites} {} threads=1 median {:.6} s max_abs {:.3e}",
            outcome.kernel, outcome.syrk.median_seconds, outcome.max_abs
        );
        let size = report
            .sizes
            .iter_mut()
            .find(|size| size.taxa == taxa && size.sites == sites)
            .expect("shape from the one-thread run");
        size.blas.push(blas_row(outcome));
    }

    let json = report::write_report(&report);
    std::fs::write(&report_path, &json).expect("write report");
    println!("wrote {}", report_path.display());
}

struct Prepared {
    taxa: usize,
    sites: usize,
    checksum: u64,
    codes: Vec<u8>,
    reference: Vec<f64>,
    reference_seconds: f64,
}

fn size_report(item: &Prepared, gpu: Vec<GpuRow>, blas: Vec<BlasRow>) -> SizeReport {
    SizeReport {
        taxa: item.taxa,
        sites: item.sites,
        checksum: item.checksum,
        preview: item.codes[..PREVIEW.len()].to_vec(),
        reference_seconds: item.reference_seconds,
        gpu,
        blas,
    }
}

fn gpu_rows(timings: &[gpu::ShaderTiming]) -> Vec<GpuRow> {
    timings
        .iter()
        .map(|timing| GpuRow {
            shader: timing.shader.clone(),
            phases: timing.phases.clone(),
            max_abs: timing.error.max_abs,
            max_rel: timing.error.max_rel,
            bit_exact: timing.error.bit_exact,
        })
        .collect()
}

fn blas_row(outcome: blas::BlasOutcome) -> BlasRow {
    BlasRow {
        kernel: outcome.kernel,
        threads: outcome.threads,
        pack: outcome.pack,
        syrk: outcome.syrk,
        max_abs: outcome.max_abs,
        max_rel: outcome.max_rel,
        bit_exact: outcome.bit_exact,
    }
}

fn write_one_thread(path: &Path) {
    eprintln!("{}", blas::library_line());
    let mut outcomes = Vec::new();
    for shape in SHAPES {
        let codes = generate_dosages(shape.taxa, shape.sites, PARK_MILLER_SEED);
        let reference = load_reference(shape.taxa, shape.sites);
        eprintln!("one-thread {}x{}", shape.taxa, shape.sites);
        for outcome in blas::time_shape(
            &BlasInput {
                taxa: shape.taxa,
                sites: shape.sites,
                codes: &codes,
                reference: &reference,
            },
            "1",
        ) {
            outcomes.push((shape.taxa, shape.sites, outcome));
        }
    }
    if let Some(parent) = path.parent() {
        std::fs::create_dir_all(parent).expect("output directory");
    }
    std::fs::write(path, blas::write_lines(&outcomes)).expect("write one-thread blas");
}

fn cache_reference(taxa: usize, sites: usize, reference: &[f64]) {
    let path = reference_path(taxa, sites);
    std::fs::create_dir_all(path.parent().unwrap()).ok();
    let mut bytes = Vec::with_capacity(reference.len() * 8);
    for value in reference {
        bytes.extend_from_slice(&value.to_ne_bytes());
    }
    std::fs::write(path, bytes).expect("cache reference");
}

fn load_reference(taxa: usize, sites: usize) -> Vec<f64> {
    let bytes = std::fs::read(reference_path(taxa, sites)).expect("cached reference");
    assert_eq!(bytes.len(), gpu_gram::pair_count(taxa) * 8);
    bytes
        .chunks_exact(8)
        .map(|chunk| f64::from_ne_bytes(chunk.try_into().unwrap()))
        .collect()
}

fn reference_path(taxa: usize, sites: usize) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join(format!("target/ref-{taxa}x{sites}.f64"))
}

fn host_report() -> HostReport {
    #[cfg(target_os = "windows")]
    {
        return windows_host();
    }
    #[cfg(not(target_os = "windows"))]
    HostReport {
        model: command("sysctl", &["-n", "hw.model"]),
        uname: command("uname", &["-srm"]),
        logical_cpus: command("sysctl", &["-n", "hw.ncpu"]),
        performance_cores: command("sysctl", &["-n", "hw.perflevel0.physicalcpu"]),
        efficiency_cores: command("sysctl", &["-n", "hw.perflevel1.physicalcpu"]),
        memory_bytes: command("sysctl", &["-n", "hw.memsize"]),
        rustc: command("rustc", &["--version"]),
    }
}

#[cfg(target_os = "windows")]
fn windows_host() -> HostReport {
    let text = command(
        "powershell.exe",
        &[
            "-NoProfile",
            "-Command",
            "[Console]::OutputEncoding = [System.Text.UTF8Encoding]::new($false); $OutputEncoding = [Console]::OutputEncoding; $cs = Get-CimInstance Win32_ComputerSystem; $os = Get-CimInstance Win32_OperatingSystem; $cpu = Get-CimInstance Win32_Processor; Write-Output $cs.Model; Write-Output ($os.Caption + ' ' + $os.Version); Write-Output $cpu.NumberOfLogicalProcessors; Write-Output $cs.TotalPhysicalMemory; Write-Output $cpu.Name",
        ],
    );
    let mut lines = text.lines().map(str::trim).filter(|line| !line.is_empty());
    let model = lines.next().unwrap_or("unknown").to_string();
    let uname = lines.next().unwrap_or("unknown").to_string();
    let logical_cpus = lines.next().unwrap_or("unknown").to_string();
    let memory_bytes = lines.next().unwrap_or("unknown").to_string();
    let cpu = lines.next().unwrap_or("unknown").to_string();
    HostReport {
        model: format!("{model}, {cpu}"),
        uname,
        logical_cpus,
        performance_cores: "unknown".to_string(),
        efficiency_cores: "unknown".to_string(),
        memory_bytes,
        rustc: command("rustc", &["--version"]),
    }
}

fn command(program: &str, args: &[&str]) -> String {
    Command::new(program)
        .args(args)
        .output()
        .ok()
        .and_then(|output| String::from_utf8(output.stdout).ok())
        .map(|text| text.trim().to_string())
        .filter(|text| !text.is_empty())
        .unwrap_or_else(|| "unknown".to_string())
}
