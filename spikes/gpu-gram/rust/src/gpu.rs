//! Headless wgpu run of `gram_f32.wgsl`. The dispatch matches `GpuGram.kt`.

use std::sync::mpsc;
use std::time::Instant;

use gpu_gram::{compare_triangles, pair_count, split_shader, PANEL_COLUMNS, ValueError};

use crate::report::Phases;

const RUNS: usize = 3;

#[repr(C)]
#[derive(Clone, Copy, bytemuck::Pod, bytemuck::Zeroable)]
struct Params {
    taxa: u32,
    columns: u32,
    stride: u32,
    pair_base: u32,
}

pub struct GpuContext {
    device: wgpu::Device,
    queue: wgpu::Queue,
    naive: wgpu::ComputePipeline,
    tiled: wgpu::ComputePipeline,
    pub adapter_name: String,
    pub vendor: u32,
    pub device_id: u32,
    pub device_type: String,
    pub backend: String,
    pub driver: String,
    pub driver_info: String,
    pub fallback: bool,
    pub limits: Vec<(String, u64)>,
    pub setup_seconds: f64,
}

pub struct ShaderTiming {
    pub shader: String,
    pub phases: Phases,
    pub error: ValueError,
}

pub fn open() -> Result<GpuContext, String> {
    let started = Instant::now();
    let instance = wgpu::Instance::new(&wgpu::InstanceDescriptor::default());
    let adapter = pollster::block_on(instance.request_adapter(&wgpu::RequestAdapterOptions {
        power_preference: wgpu::PowerPreference::HighPerformance,
        force_fallback_adapter: false,
        compatible_surface: None,
    }))
    .map_err(|error| format!("request_adapter: {error:?}"))?;
    let info = adapter.get_info();
    let available = adapter.limits();
    let (device, queue) = pollster::block_on(adapter.request_device(&wgpu::DeviceDescriptor::default()))
        .map_err(|error| format!("request_device: {error:?}"))?;
    device.on_uncaptured_error(Box::new(|error| panic!("wgpu: {error:?}")));
    let (naive_source, tiled_source) = {
        let path = gpu_gram::shader_path();
        let source = std::fs::read_to_string(&path).map_err(|error| format!("reading {}: {error}", path.display()))?;
        split_shader(&source)
    };
    let naive_module = device.create_shader_module(wgpu::ShaderModuleDescriptor {
        label: Some("gram_naive"),
        source: wgpu::ShaderSource::Wgsl(naive_source.into()),
    });
    let tiled_module = device.create_shader_module(wgpu::ShaderModuleDescriptor {
        label: Some("gram_tiled"),
        source: wgpu::ShaderSource::Wgsl(tiled_source.into()),
    });
    let naive = device.create_compute_pipeline(&wgpu::ComputePipelineDescriptor {
        label: Some("gram_naive"),
        layout: None,
        module: &naive_module,
        entry_point: Some("gram_naive"),
        compilation_options: Default::default(),
        cache: None,
    });
    let tiled = device.create_compute_pipeline(&wgpu::ComputePipelineDescriptor {
        label: Some("gram_tiled"),
        layout: None,
        module: &tiled_module,
        entry_point: Some("gram_tiled"),
        compilation_options: Default::default(),
        cache: None,
    });
    let granted = device.limits();
    let fallback = matches!(info.device_type, wgpu::DeviceType::Cpu | wgpu::DeviceType::VirtualGpu);
    let limits = vec![
        ("adapter.maxBufferSize".to_string(), available.max_buffer_size),
        ("adapter.maxStorageBufferBindingSize".to_string(), available.max_storage_buffer_binding_size as u64),
        ("adapter.maxComputeWorkgroupSizeX".to_string(), available.max_compute_workgroup_size_x as u64),
        ("adapter.maxComputeInvocationsPerWorkgroup".to_string(), available.max_compute_invocations_per_workgroup as u64),
        ("adapter.maxComputeWorkgroupsPerDimension".to_string(), available.max_compute_workgroups_per_dimension as u64),
        ("device.maxBufferSize".to_string(), granted.max_buffer_size),
        ("device.maxStorageBufferBindingSize".to_string(), granted.max_storage_buffer_binding_size as u64),
        ("device.maxComputeWorkgroupSizeX".to_string(), granted.max_compute_workgroup_size_x as u64),
        ("device.maxComputeInvocationsPerWorkgroup".to_string(), granted.max_compute_invocations_per_workgroup as u64),
        ("device.maxComputeWorkgroupsPerDimension".to_string(), granted.max_compute_workgroups_per_dimension as u64),
    ];
    Ok(GpuContext {
        device,
        queue,
        naive,
        tiled,
        adapter_name: info.name,
        vendor: info.vendor,
        device_id: info.device,
        device_type: format!("{:?}", info.device_type),
        backend: format!("{:?}", info.backend),
        driver: info.driver,
        driver_info: info.driver_info,
        fallback,
        limits,
        setup_seconds: elapsed(started),
    })
}

pub fn time_shape(gpu: &GpuContext, codes: &[u8], taxa: usize, sites: usize, reference: &[f64]) -> Result<[ShaderTiming; 2], String> {
    let pairs = pair_count(taxa);
    let triangle_bytes = (pairs * std::mem::size_of::<f32>()) as u64;
    let read_bytes = align_up(triangle_bytes, wgpu::MAP_ALIGNMENT);
    let panel_bytes = (taxa * PANEL_COLUMNS * std::mem::size_of::<f32>()) as u64;
    let panel = gpu.device.create_buffer(&wgpu::BufferDescriptor {
        label: Some("panel"),
        size: panel_bytes,
        usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_DST,
        mapped_at_creation: false,
    });
    let triangle = gpu.device.create_buffer(&wgpu::BufferDescriptor {
        label: Some("triangle"),
        size: triangle_bytes,
        usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_SRC | wgpu::BufferUsages::COPY_DST,
        mapped_at_creation: false,
    });
    let read = gpu.device.create_buffer(&wgpu::BufferDescriptor {
        label: Some("readback"),
        size: read_bytes,
        usage: wgpu::BufferUsages::MAP_READ | wgpu::BufferUsages::COPY_DST,
        mapped_at_creation: false,
    });
    let uniform = gpu.device.create_buffer(&wgpu::BufferDescriptor {
        label: Some("params"),
        size: std::mem::size_of::<Params>() as u64,
        usage: wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST,
        mapped_at_creation: false,
    });
    let naive_group = bind_group(&gpu.device, &gpu.naive, &panel, &triangle, &uniform);
    let tiled_group = bind_group(&gpu.device, &gpu.tiled, &panel, &triangle, &uniform);
    let naive_spans = dispatch_spans(pairs as u32, 256);
    let tiled_spans = dispatch_spans(upper_tiles(taxa as u32), 1);
    let naive = time_shader(gpu, &gpu.naive, &naive_group, &panel, &triangle, &read, &uniform, codes, taxa, sites, triangle_bytes, &naive_spans, reference)?;
    let tiled = time_shader(gpu, &gpu.tiled, &tiled_group, &panel, &triangle, &read, &uniform, codes, taxa, sites, triangle_bytes, &tiled_spans, reference)?;
    Ok([
        ShaderTiming { shader: "naive".to_string(), phases: naive.0, error: naive.1 },
        ShaderTiming { shader: "tiled".to_string(), phases: tiled.0, error: tiled.1 },
    ])
}

fn bind_group(
    device: &wgpu::Device,
    pipeline: &wgpu::ComputePipeline,
    panel: &wgpu::Buffer,
    triangle: &wgpu::Buffer,
    uniform: &wgpu::Buffer,
) -> wgpu::BindGroup {
    device.create_bind_group(&wgpu::BindGroupDescriptor {
        label: None,
        layout: &pipeline.get_bind_group_layout(0),
        entries: &[
            wgpu::BindGroupEntry { binding: 0, resource: panel.as_entire_binding() },
            wgpu::BindGroupEntry { binding: 1, resource: triangle.as_entire_binding() },
            wgpu::BindGroupEntry { binding: 2, resource: uniform.as_entire_binding() },
        ],
    })
}

fn time_shader(
    gpu: &GpuContext,
    pipeline: &wgpu::ComputePipeline,
    bind_group: &wgpu::BindGroup,
    panel_buffer: &wgpu::Buffer,
    triangle_buffer: &wgpu::Buffer,
    read_buffer: &wgpu::Buffer,
    uniform: &wgpu::Buffer,
    codes: &[u8],
    taxa: usize,
    sites: usize,
    triangle_bytes: u64,
    spans: &[(u32, u32)],
    reference: &[f64],
) -> Result<(Phases, ValueError), String> {
    let mut staging = Vec::with_capacity(RUNS);
    let mut write_buffer = Vec::with_capacity(RUNS);
    let mut dispatch = Vec::with_capacity(RUNS);
    let mut map_async = Vec::with_capacity(RUNS);
    let mut copy_back = Vec::with_capacity(RUNS);
    let mut cpu_accumulation = Vec::with_capacity(RUNS);
    let mut end_to_end = Vec::with_capacity(RUNS);
    let mut last = Vec::new();
    let mut panel = vec![0.0f32; taxa * PANEL_COLUMNS];
    for iteration in 0..=RUNS {
        let mut staging_seconds = 0.0;
        let mut write_seconds = 0.0;
        let mut dispatch_seconds = 0.0;
        let mut wrote_any = false;
        let mut site = 0;
        while site < sites {
            let columns = PANEL_COLUMNS.min(sites - site);
            let stage_mark = Instant::now();
            fill_panel(codes, taxa, site, columns, &mut panel);
            staging_seconds += elapsed(stage_mark);
            let write_mark = Instant::now();
            let bytes = taxa * columns * std::mem::size_of::<f32>();
            gpu.queue.write_buffer(panel_buffer, 0, bytemuck::cast_slice(&panel[..taxa * columns]));
            debug_assert_eq!(bytes % 4, 0);
            write_seconds += elapsed(write_mark);
            for (span_index, &(pair_base, groups)) in spans.iter().enumerate() {
                let clear = !wrote_any && span_index == 0;
                let dispatch_mark = Instant::now();
                let params = Params {
                    taxa: taxa as u32,
                    columns: columns as u32,
                    stride: columns as u32,
                    pair_base,
                };
                gpu.queue.write_buffer(uniform, 0, bytemuck::bytes_of(&params));
                let mut encoder = gpu.device.create_command_encoder(&wgpu::CommandEncoderDescriptor { label: None });
                if clear {
                    encoder.clear_buffer(triangle_buffer, 0, None);
                }
                {
                    let mut pass = encoder.begin_compute_pass(&wgpu::ComputePassDescriptor {
                        label: None,
                        timestamp_writes: None,
                    });
                    pass.set_pipeline(pipeline);
                    pass.set_bind_group(0, bind_group, &[]);
                    pass.dispatch_workgroups(groups, 1, 1);
                }
                let submission = gpu.queue.submit([encoder.finish()]);
                gpu.device
                    .poll(wgpu::PollType::wait_for(submission))
                    .map_err(|error| format!("dispatch poll: {error:?}"))?;
                dispatch_seconds += elapsed(dispatch_mark);
                wrote_any = true;
            }
            site += columns;
        }
        let map_mark = Instant::now();
        let mut encoder = gpu.device.create_command_encoder(&wgpu::CommandEncoderDescriptor { label: None });
        encoder.copy_buffer_to_buffer(triangle_buffer, 0, read_buffer, 0, triangle_bytes);
        let submission = gpu.queue.submit([encoder.finish()]);
        gpu.device
            .poll(wgpu::PollType::wait_for(submission))
            .map_err(|error| format!("readback poll: {error:?}"))?;
        let (sender, receiver) = mpsc::channel();
        read_buffer.slice(..).map_async(wgpu::MapMode::Read, move |result| {
            let _ = sender.send(result);
        });
        gpu.device
            .poll(wgpu::PollType::wait())
            .map_err(|error| format!("map poll: {error:?}"))?;
        receiver
            .recv()
            .map_err(|_| "map_async callback was dropped".to_string())?
            .map_err(|error| format!("map_async: {error:?}"))?;
        let map_seconds = elapsed(map_mark);
        let copy_mark = Instant::now();
        let mapped = read_buffer.slice(..).get_mapped_range();
        let floats: &[f32] = bytemuck::cast_slice(&mapped[..triangle_bytes as usize]);
        let triangle = floats.to_vec();
        drop(mapped);
        read_buffer.unmap();
        let copy_seconds = elapsed(copy_mark);
        if iteration == 0 {
            continue;
        }
        let cpu_seconds = 0.0;
        staging.push(staging_seconds);
        write_buffer.push(write_seconds);
        dispatch.push(dispatch_seconds);
        map_async.push(map_seconds);
        copy_back.push(copy_seconds);
        cpu_accumulation.push(cpu_seconds);
        end_to_end.push(staging_seconds + write_seconds + dispatch_seconds + map_seconds + copy_seconds + cpu_seconds);
        last = triangle;
    }
    let error = compare_f32(&last, reference);
    Ok((
        Phases {
            staging: crate::report::stats_of(&staging),
            write_buffer: crate::report::stats_of(&write_buffer),
            dispatch: crate::report::stats_of(&dispatch),
            map_async: crate::report::stats_of(&map_async),
            copy_back: crate::report::stats_of(&copy_back),
            cpu_accumulation: crate::report::stats_of(&cpu_accumulation),
            end_to_end: crate::report::stats_of(&end_to_end),
        },
        error,
    ))
}

fn compare_f32(actual: &[f32], reference: &[f64]) -> ValueError {
    assert_eq!(actual.len(), reference.len(), "triangles differ");
    let promoted: Vec<f64> = actual.iter().copied().map(f64::from).collect();
    compare_triangles(&promoted, reference)
}

fn fill_panel(codes: &[u8], taxa: usize, site_start: usize, columns: usize, panel: &mut [f32]) {
    for taxon in 0..taxa {
        let row = taxon * columns;
        for column in 0..columns {
            let code = codes[(site_start + column) * taxa + taxon];
            panel[row + column] = if code == gpu_gram::MISSING_CODE { 0.0 } else { code as f32 };
        }
    }
}

fn dispatch_spans(item_count: u32, items_per_group: u32) -> Vec<(u32, u32)> {
    let mut spans = Vec::new();
    let mut base = 0u32;
    while base < item_count {
        let groups = ((item_count - base + items_per_group - 1) / items_per_group).min(65_535);
        spans.push((base, groups));
        base += groups * items_per_group;
    }
    spans
}

fn upper_tiles(taxa: u32) -> u32 {
    let tiles = (taxa + 15) / 16;
    tiles * (tiles + 1) / 2
}

fn align_up(value: u64, align: u64) -> u64 {
    (value + align - 1) & !(align - 1)
}

fn elapsed(started: Instant) -> f64 {
    started.elapsed().as_secs_f64()
}
