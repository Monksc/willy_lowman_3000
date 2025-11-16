use willy_loman_3000::*;

pub fn filter(mut points: Vec<Point>, n: usize) -> Vec<Point> {
    points.truncate(n);
    points
}

pub fn berlin() -> Vec<Point> {
    vec![
        (565.0, 575.0).into(),
        (25.0, 185.0).into(),
        (345.0, 750.0).into(),
        (945.0, 685.0).into(),
        (845.0, 655.0).into(),
        (880.0, 660.0).into(),
        (25.0, 230.0).into(),
        (525.0, 1000.0).into(),
        (580.0, 1175.0).into(),
        (650.0, 1130.0).into(),
        (1605.0, 620.0).into(),
        (1220.0, 580.0).into(),
        (1465.0, 200.0).into(),
        (1530.0, 5.0).into(),
        (845.0, 680.0).into(),
        (725.0, 370.0).into(),
        (145.0, 665.0).into(),
        (415.0, 635.0).into(),
        (510.0, 875.0).into(),
        (560.0, 365.0).into(),
        (300.0, 465.0).into(),
        (520.0, 585.0).into(),
        (480.0, 415.0).into(),
        (835.0, 625.0).into(),
        (975.0, 580.0).into(),
        (1215.0, 245.0).into(),
        (1320.0, 315.0).into(),
        (1250.0, 400.0).into(),
        (660.0, 180.0).into(),
        (410.0, 250.0).into(),
        (420.0, 555.0).into(),
        (575.0, 665.0).into(),
        (1150.0, 1160.0).into(),
        (700.0, 580.0).into(),
        (685.0, 595.0).into(),
        (685.0, 610.0).into(),
        (770.0, 610.0).into(),
        (795.0, 645.0).into(),
        (720.0, 635.0).into(),
        (760.0, 650.0).into(),
        (475.0, 960.0).into(),
        (95.0, 260.0).into(),
        (875.0, 920.0).into(),
        (700.0, 500.0).into(),
        (555.0, 815.0).into(),
        (830.0, 485.0).into(),
        (1170.0, 65.0).into(),
        (830.0, 610.0).into(),
        (605.0, 625.0).into(),
        (595.0, 360.0).into(),
        (1340.0, 725.0).into(),
        (1740.0, 245.0).into(),
    ]
}

pub fn willie_loman_held_karp(
    restarts: usize,
    max_k: usize,
    y_partition: f32,
) -> Box<dyn Fn(Vec<Point>) -> f32> {
    Box::from(move |points: Vec<Point>| -> f32 {
        let mut tsp: TSP = make_graph_points(&points).into();

        tsp.update_solution(snake_pattern(&points, y_partition));

        tsp.chained_lin_kernighan(restarts, max_k);
        tsp.held_karp_solution();

        tsp.update_cost();
        tsp.cost()
    })
}

pub fn willie_loman(
    restarts: usize,
    max_k: usize,
    y_partition: f32,
) -> Box<dyn Fn(Vec<Point>) -> f32> {
    Box::from(move |points: Vec<Point>| -> f32 {
        let mut tsp: TSP = make_graph_points(&points).into();

        tsp.update_solution(snake_pattern(&points, y_partition));

        tsp.chained_lin_kernighan(restarts, max_k);
        tsp.update_cost();
        tsp.cost()
    })
}

pub fn concorde() -> Box<dyn Fn(Vec<Point>) -> f32> {
    use concorde_rs::*;

    struct Node(i32, i32);

    impl Distance for Node {
        fn calc_shortest_dist(&self, other: &Self) -> u32 {
            self.0.abs_diff(other.0) + self.1.abs_diff(other.1)
        }
    }

    Box::from(move |points: Vec<Point>| -> f32 {
        let tsp: TSP = make_graph_points(&points).into();

        let nodes: Vec<Node> = points
            .into_iter()
            .map(|p| Node(p.x as i32, p.y as i32))
            .collect();
        let dist_mat = LowerDistanceMatrix::from(nodes.as_ref());
        let solution = solver::tsp_hk(&dist_mat).unwrap();

        let solution: Vec<usize> = solution
            .tour
            .into_iter()
            .map(|index| index as usize)
            .collect();
        tsp.cost_of_solution(&solution)
    })
}

macro_rules! bench_tsp {
    ($name:ident, $data:expr, $tspalgo:expr) => {
        fn $name() {
            eprintln!("Name: {} Cost: {}", stringify!($name), $tspalgo($data));
        }
    };
}
pub(crate) use bench_tsp;
