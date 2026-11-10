//! JSON for the harness report. Written by hand so the crate stays on wgpu, pollster, and bytemuck.

use std::fmt::Write as _;

#[derive(Clone)]
pub struct Stats {
    pub min_seconds: f64,
    pub median_seconds: f64,
    pub samples: Vec<f64>,
}

pub fn stats_of(samples: &[f64]) -> Stats {
    assert!(!samples.is_empty(), "no samples");
    let mut sorted = samples.to_vec();
    sorted.sort_by(|left, right| left.partial_cmp(right).unwrap());
    let median_at = (sorted.len() + 1) / 2 - 1;
    Stats {
        min_seconds: sorted[0],
        median_seconds: sorted[median_at],
        samples: samples.to_vec(),
    }
}

#[derive(Clone)]
pub struct Phases {
    pub staging: Stats,
    pub write_buffer: Stats,
    pub dispatch: Stats,
    pub map_async: Stats,
    pub copy_back: Stats,
    pub cpu_accumulation: Stats,
    pub end_to_end: Stats,
}

pub struct GpuRow {
    pub shader: String,
    pub phases: Phases,
    pub max_abs: f64,
    pub max_rel: f64,
    pub bit_exact: bool,
}

pub struct BlasRow {
    pub kernel: String,
    pub threads: String,
    pub pack: Stats,
    pub syrk: Stats,
    pub max_abs: f64,
    pub max_rel: f64,
    pub bit_exact: bool,
}

pub struct SizeReport {
    pub taxa: usize,
    pub sites: usize,
    pub checksum: u64,
    pub preview: Vec<u8>,
    pub reference_seconds: f64,
    pub gpu: Vec<GpuRow>,
    pub blas: Vec<BlasRow>,
}

pub struct HostReport {
    pub model: String,
    pub uname: String,
    pub logical_cpus: String,
    pub performance_cores: String,
    pub efficiency_cores: String,
    pub memory_bytes: String,
    pub rustc: String,
}

pub struct AdapterReport {
    pub name: String,
    pub vendor: u32,
    pub device: u32,
    pub device_type: String,
    pub backend: String,
    pub driver: String,
    pub driver_info: String,
    pub fallback: bool,
    pub limits: Vec<(String, u64)>,
}

pub struct Report {
    pub host: HostReport,
    pub adapter: Option<AdapterReport>,
    pub gpu_error: Option<String>,
    pub setup_seconds: Option<f64>,
    pub sizes: Vec<SizeReport>,
}

pub fn write_report(report: &Report) -> String {
    let mut out = String::new();
    out.push_str("{\n");
    out.push_str("  \"host\": ");
    write_host(&mut out, &report.host);
    out.push_str(",\n  \"adapter\": ");
    match &report.adapter {
        Some(adapter) => write_adapter(&mut out, adapter),
        None => out.push_str("null"),
    }
    out.push_str(",\n  \"gpuError\": ");
    write_opt_string(&mut out, report.gpu_error.as_deref());
    out.push_str(",\n  \"setupSeconds\": ");
    match report.setup_seconds {
        Some(seconds) => out.push_str(&json_seconds(seconds)),
        None => out.push_str("null"),
    }
    out.push_str(",\n  \"sizes\": [\n");
    for (index, size) in report.sizes.iter().enumerate() {
        if index > 0 {
            out.push_str(",\n");
        }
        write_size(&mut out, size);
    }
    out.push_str("\n  ]\n}\n");
    out
}

fn write_host(out: &mut String, host: &HostReport) {
    let _ = write!(
        out,
        "{{\n    \"model\": {},\n    \"uname\": {},\n    \"logicalCpus\": {},\n    \"performanceCores\": {},\n    \"efficiencyCores\": {},\n    \"memoryBytes\": {},\n    \"rustc\": {}\n  }}",
        json_string(&host.model),
        json_string(&host.uname),
        json_string(&host.logical_cpus),
        json_string(&host.performance_cores),
        json_string(&host.efficiency_cores),
        json_string(&host.memory_bytes),
        json_string(&host.rustc),
    );
}

fn write_adapter(out: &mut String, adapter: &AdapterReport) {
    let _ = write!(
        out,
        "{{\n    \"name\": {},\n    \"vendor\": {},\n    \"device\": {},\n    \"deviceType\": {},\n    \"backend\": {},\n    \"driver\": {},\n    \"driverInfo\": {},\n    \"fallback\": {},\n    \"limits\": {{",
        json_string(&adapter.name),
        adapter.vendor,
        adapter.device,
        json_string(&adapter.device_type),
        json_string(&adapter.backend),
        json_string(&adapter.driver),
        json_string(&adapter.driver_info),
        adapter.fallback,
    );
    for (index, (name, value)) in adapter.limits.iter().enumerate() {
        if index > 0 {
            out.push(',');
        }
        let _ = write!(out, "\n      {}: {value}", json_string(name));
    }
    out.push_str("\n    }\n  }");
}

fn write_size(out: &mut String, size: &SizeReport) {
    let _ = write!(
        out,
        "    {{\n      \"taxa\": {},\n      \"sites\": {},\n      \"checksum\": \"{}\",\n      \"preview\": [",
        size.taxa, size.sites, size.checksum
    );
    for (index, code) in size.preview.iter().enumerate() {
        if index > 0 {
            out.push_str(", ");
        }
        let _ = write!(out, "{code}");
    }
    let _ = write!(
        out,
        "],\n      \"referenceSeconds\": {},\n      \"gpu\": [",
        json_seconds(size.reference_seconds)
    );
    for (index, row) in size.gpu.iter().enumerate() {
        if index > 0 {
            out.push(',');
        }
        write_gpu(out, row);
    }
    out.push_str("\n      ],\n      \"blas\": [");
    for (index, row) in size.blas.iter().enumerate() {
        if index > 0 {
            out.push(',');
        }
        write_blas(out, row);
    }
    out.push_str("\n      ]\n    }");
}

fn write_gpu(out: &mut String, row: &GpuRow) {
    let _ = write!(
        out,
        "\n        {{\n          \"shader\": {},\n          \"transfer\": \"write_buffer\",\n          \"accumulation\": \"f32\",\n          \"panel\": \"dosage\",\n          \"phases\": ",
        json_string(&row.shader)
    );
    write_phases(out, &row.phases);
    let _ = write!(
        out,
        ",\n          \"maxAbs\": {},\n          \"maxRel\": {},\n          \"dosageBitExact\": {}\n        }}",
        json_error(row.max_abs),
        json_error(row.max_rel),
        row.bit_exact,
    );
}

fn write_blas(out: &mut String, row: &BlasRow) {
    let _ = write!(
        out,
        "\n        {{\n          \"kernel\": {},\n          \"threads\": {},\n          \"pack\": ",
        json_string(&row.kernel),
        json_string(&row.threads),
    );
    write_stats(out, &row.pack);
    out.push_str(",\n          \"syrk\": ");
    write_stats(out, &row.syrk);
    let _ = write!(
        out,
        ",\n          \"maxAbs\": {},\n          \"maxRel\": {},\n          \"bitExact\": {}\n        }}",
        json_error(row.max_abs),
        json_error(row.max_rel),
        row.bit_exact,
    );
}

fn write_phases(out: &mut String, phases: &Phases) {
    out.push('{');
    for (index, (name, stats)) in [
        ("staging", &phases.staging),
        ("writeBuffer", &phases.write_buffer),
        ("dispatch", &phases.dispatch),
        ("mapAsync", &phases.map_async),
        ("copyBack", &phases.copy_back),
        ("cpuAccumulation", &phases.cpu_accumulation),
        ("endToEnd", &phases.end_to_end),
    ]
    .into_iter()
    .enumerate()
    {
        if index > 0 {
            out.push(',');
        }
        let _ = write!(out, "\n            {}: ", json_string(name));
        write_stats(out, stats);
    }
    out.push_str("\n          }");
}

fn write_stats(out: &mut String, stats: &Stats) {
    let _ = write!(
        out,
        "{{\"minSeconds\": {}, \"medianSeconds\": {}, \"samplesSeconds\": [",
        json_seconds(stats.min_seconds),
        json_seconds(stats.median_seconds),
    );
    for (index, sample) in stats.samples.iter().enumerate() {
        if index > 0 {
            out.push_str(", ");
        }
        out.push_str(&json_seconds(*sample));
    }
    out.push_str("]}");
}

fn write_opt_string(out: &mut String, value: Option<&str>) {
    match value {
        Some(text) => out.push_str(&json_string(text)),
        None => out.push_str("null"),
    }
}

fn json_string(value: &str) -> String {
    let mut out = String::from("\"");
    for ch in value.chars() {
        match ch {
            '\\' => out.push_str("\\\\"),
            '"' => out.push_str("\\\""),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            other => out.push(other),
        }
    }
    out.push('"');
    out
}

fn json_seconds(value: f64) -> String {
    if value == 0.0 {
        "0.0".to_string()
    } else {
        format!("{value:.9}")
    }
}

fn json_error(value: f64) -> String {
    if value == 0.0 {
        "0.0".to_string()
    } else {
        format!("{value:.16e}")
    }
}
