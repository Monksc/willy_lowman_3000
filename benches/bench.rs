#[path = "data_algo.rs"]
mod data_algo;
use data_algo::*;

data_algo::bench_tsp!(berlin_willie_loman, berlin(), willie_loman(256, 1, 250.0));
data_algo::bench_tsp!(
    berlin_willie_loman_held_karp,
    filter(berlin(), 22),
    willie_loman_held_karp(256, 1, 250.0)
);
data_algo::bench_tsp!(berlin_concorde_rs, berlin(), concorde());

iai::main!(
    berlin_willie_loman,
    berlin_concorde_rs,
    berlin_willie_loman_held_karp
);
