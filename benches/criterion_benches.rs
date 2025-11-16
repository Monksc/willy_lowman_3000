#[path = "data_algo.rs"]
mod data_algo;
use data_algo::*;

use std::time::Duration;

use criterion::{criterion_group, criterion_main, Criterion};

data_algo::bench_tsp!(
    berlin_willie_loman,
    berlin(),
    willie_loman(256, 1000, 250.0)
);

data_algo::bench_tsp!(
    berlin_willie_loman_held_karp,
    filter(berlin(), 20),
    willie_loman_held_karp(256, 1, 250.0)
);

data_algo::bench_tsp!(berlin_concorde_rs, berlin(), concorde());

fn criterion_benchmark(c: &mut Criterion) {
    c.bench_function("berlin_willie_loman", |b| b.iter(|| berlin_willie_loman()));
    c.bench_function("berlin_willie_loman_held_karp", |b| {
        b.iter(|| berlin_willie_loman_held_karp())
    });
    c.bench_function("berlin_concorde_rs", |b| b.iter(|| berlin_concorde_rs()));
}

criterion_group! {
    name = benches;
    config = Criterion::default().sample_size(10).measurement_time(Duration::from_secs(100));
    targets = criterion_benchmark
}
criterion_main!(benches);
