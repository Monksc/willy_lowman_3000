#[derive(Debug, Clone)]
pub struct TSP {
    weights: Vec<Vec<f32>>,
    _closest_index: Vec<Vec<usize>>,
    solution: Vec<usize>,
    vertex_to_solution_index: Vec<usize>,
    cost: f32,
}

impl From<Vec<Vec<f64>>> for TSP {
    fn from(value: Vec<Vec<f64>>) -> Self {
        value
            .into_iter()
            .map(|x| x.into_iter().map(|x| x as f32).collect())
            .collect::<Vec<Vec<f32>>>()
            .into()
    }
}

impl From<Vec<Vec<f32>>> for TSP {
    fn from(value: Vec<Vec<f32>>) -> Self {
        let solution: Vec<usize> = (0..value.len()).collect();
        let closest_index = (0..value.len())
            .map(|index| {
                let mut row: Vec<(f32, usize)> = value[index]
                    .iter()
                    .enumerate()
                    .map(|(inner_index, value)| (*value, inner_index))
                    .collect();

                row.sort_by(|(a, _), (b, _)| a.partial_cmp(b).unwrap());

                row.into_iter().map(|(_, index)| index).collect()
            })
            .collect();
        let mut s = Self {
            weights: value,
            _closest_index: closest_index,
            solution: Vec::new(),
            vertex_to_solution_index: solution.clone(),
            cost: 0.0,
        };

        s.cost = s.cost_of_solution(&solution);
        s.solution = solution;
        s
    }
}

impl TSP {
    pub fn cost_of_solution(&self, solution: &[usize]) -> f32 {
        let mut c = 0.0;
        for i in 1..solution.len() {
            c += self.weights[solution[i - 1]][solution[i]];
        }
        c
    }

    pub fn cost(&self) -> f32 {
        self.cost
    }

    pub fn copy_reset_solution(&self) -> Self {
        let mut weights = self.weights.clone();
        for i in 0..self.solution.len() {
            for j in 0..self.solution.len() {
                weights[i][j] = self.weights[self.solution[i]][self.solution[j]];
            }
        }

        weights.into()
    }

    pub fn update_solution(&mut self, solution: Vec<usize>) {
        let c = self.cost_of_solution(&solution);
        if c < self.cost {
            self.cost = c;
            self.solution = solution;
            for i in 0..self.solution.len() {
                self.vertex_to_solution_index[self.solution[i]] = i;
            }
        }
    }

    pub fn update_cost(&mut self) {
        self.cost = self.cost_of_solution(&self.solution);
    }

    pub fn cycle_cost(&self) -> f32 {
        self.cost + self.weights[self.solution[self.solution.len() - 1]][self.solution[0]]
    }

    pub fn solution(&self) -> &Vec<usize> {
        &self.solution
    }

    #[inline]
    fn prev_index(&self, solution_index: usize) -> usize {
        // self.solution[(solution_index + self.weights.len() - 1) % self.weights.len()]
        self.solution[solution_index - 1]
    }

    #[inline]
    fn next_index(&self, solution_index: usize) -> usize {
        // self.solution[(solution_index + 1) % self.weights.len()]
        self.solution[solution_index + 1]
    }

    fn cost_index(&self, solution_index: usize) -> f32 {
        let mut cost = 0.0;
        if solution_index > 0 {
            cost += self.weights[self.prev_index(solution_index)][self.solution[solution_index]];
        }

        if solution_index + 1 < self.weights.len() {
            cost += self.weights[self.solution[solution_index]][self.next_index(solution_index)]
        }

        cost
    }

    fn swap(&mut self, j_solution_index: usize, k_solution_index: usize) {
        let pre_cost = self.cost_index(j_solution_index) + self.cost_index(k_solution_index);
        let j = self.solution[j_solution_index];
        let k = self.solution[k_solution_index];
        self.solution.swap(j_solution_index, k_solution_index);
        self.vertex_to_solution_index.swap(j, k);
        let post_cost = self.cost_index(j_solution_index) + self.cost_index(k_solution_index);
        self.cost += post_cost - pre_cost;
        debug_assert!((self.cost - self.cost_of_solution(&self.solution)).abs() < 0.1);
    }

    fn try_swap(&mut self, j_solution_index: usize, k_solution_index: usize) -> bool {
        let pre_cost = self.cost;
        self.swap(j_solution_index, k_solution_index);
        let post_cost = self.cost;

        if pre_cost <= post_cost {
            self.swap(j_solution_index, k_solution_index);
            false
        } else {
            true
        }
    }

    pub fn random_swap(&mut self) {
        let j = rand::random_range(0..self.weights.len() - 1);
        let k = rand::random_range(j + 1..self.weights.len());
        self.swap(j, k);
    }

    pub fn random_swap_n(&mut self, n: usize) {
        for _ in 0..n {
            self.random_swap();
        }
    }

    pub fn chained_lin_kernighan(&mut self, restarts: usize, max_k: usize) {
        if self.weights.len() < 2 {
            return;
        }
        self.lin_kernighan(max_k);

        use rand::Rng;

        let mut rng = rand::rng();
        for _ in 0..restarts {
            let mut new_solution = self.clone();
            let random_numb: f32 = rng.random();
            let random_numb = random_numb * random_numb;
            let random_numb = (random_numb * self.weights.len() as f32).ceil() as usize;

            new_solution.random_swap_n(if random_numb < 2 { 2 } else { random_numb });
            new_solution.lin_kernighan(max_k);

            if new_solution.cost < self.cost {
                *self = new_solution;
            }
        }
    }

    pub fn lin_kernighan(&mut self, _max_k: usize) {
        'outerloop: loop {
            for index in 0..self.weights.len() {
                for second_index in (index + 1)..self.weights.len() {
                    if self.try_swap(
                        index,
                        second_index, // self.vertex_to_solution_index[self.closest_index[index][second_index]],
                    ) {
                        continue 'outerloop;
                    }
                }
            }

            break;
        }
    }

    pub fn held_karp_solution(&mut self) {
        let mut s = self.copy_reset_solution();
        s.held_karp();

        let solution = s
            .solution
            .into_iter()
            .map(|index| self.vertex_to_solution_index[index])
            .collect();

        self.update_solution(solution);
    }

    pub fn held_karp(&mut self) {
        let mut memo = std::collections::HashMap::new();
        self.dp(1, 0, &mut memo);
        let solution = self.reconstruct_path(&mut memo);
        self.update_solution(solution);
    }

    fn reconstruct_path(
        &self,
        memo: &mut std::collections::HashMap<(usize, usize), (f32, usize)>,
    ) -> Vec<usize> {
        let mut path = Vec::with_capacity(self.weights.len());
        let mut visited = 1;
        let mut last = 0;

        while path.len() < self.weights.len() {
            path.push(last);
            if let Some(&(_, next)) = memo.get(&(visited, last)) {
                visited |= 1 << next;
                last = next;
            } else {
                break;
            }
        }

        path.push(0); // Return to the start city
        path
    }

    fn dp(
        &mut self,
        visited: usize,
        last: usize,
        memo: &mut std::collections::HashMap<(usize, usize), (f32, usize)>,
    ) -> f32 {
        // If all cities are visited, return cost to go back to the start
        if visited == (1 << self.weights.len()) - 1 {
            return self.weights[last][0];
        }

        // Check if already computed
        if let Some(&(cost, _)) = memo.get(&(visited, last)) {
            return cost;
        }

        let mut min_cost = f32::INFINITY;
        let mut best_prev = usize::MAX;

        for next in 0..self.weights.len() {
            if visited & (1 << next) == 0 {
                // Try visiting 'next' city
                let new_visited = visited | (1 << next);
                let new_cost = self.weights[last][next] + self.dp(new_visited, next, memo);
                if new_cost < min_cost {
                    min_cost = new_cost;
                    best_prev = next;
                }
            }
        }

        // Memoize and return result
        memo.insert((visited, last), (min_cost, best_prev));
        min_cost
    }
}

#[derive(Debug, Clone)]
pub struct Point {
    pub x: f32,
    pub y: f32,
}

impl From<(f64, f64)> for Point {
    fn from(value: (f64, f64)) -> Self {
        Self {
            x: value.0 as f32,
            y: value.1 as f32,
        }
    }
}

impl From<(f32, f32)> for Point {
    fn from(value: (f32, f32)) -> Self {
        Self {
            x: value.0,
            y: value.1,
        }
    }
}

impl Point {
    pub fn distance(&self, other: &Point) -> f32 {
        let difx = self.x - other.x;
        let dify = self.y - other.y;
        (difx * difx + dify * dify).sqrt()
    }
}

pub fn snake_pattern(points: &[Point], y_partition: f32) -> Vec<usize> {
    let mut points = points
        .iter()
        .map(|point| point.clone())
        .enumerate()
        .collect::<Vec<(usize, Point)>>();

    points.sort_by(|(_, l), (_, r)| {
        let ly_level = (l.y / y_partition).floor() as usize;
        let ry_level = (r.y / y_partition).floor() as usize;

        if ly_level == ry_level {
            return if (ly_level % 2 == 0) == (l.y < r.y) {
                std::cmp::Ordering::Less
            } else {
                std::cmp::Ordering::Greater
            };
        }

        if ly_level < ry_level {
            std::cmp::Ordering::Less
        } else {
            std::cmp::Ordering::Greater
        }
    });

    points.into_iter().map(|(index, _)| index).collect()
}

pub fn make_graph_points(points: &Vec<Point>) -> Vec<Vec<f32>> {
    let mut graph: Vec<Vec<f32>> = (0..points.len())
        .map(|_| (0..points.len()).map(|_| 0.0).collect())
        .collect();

    for i in 0..points.len() {
        for j in (i + 1)..points.len() {
            let d = points[i].distance(&points[j]);
            graph[i][j] = d;
            graph[j][i] = d;
        }
    }

    graph
}

pub fn make_graph<T: Into<Point> + Clone>(points: &Vec<T>) -> Vec<Vec<f32>> {
    make_graph_points(
        &points
            .iter()
            .map(|point| point.clone().into())
            .collect::<Vec<Point>>(),
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    fn berlin() -> Vec<Point> {
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

    #[test]
    fn simple_test_1() {
        let mut tsp1: TSP = vec![
            vec![0.0, 1.0, 10.0],
            vec![1.0, 0.0, 10.0],
            vec![10.0, 10.0, 0.0],
        ]
        .into();

        assert_eq!(tsp1.cost, 11.0);

        tsp1.chained_lin_kernighan(1_000, 3);

        assert_eq!(tsp1.cost, 11.0);
        assert_eq!(tsp1.solution, vec![0, 1, 2]);
    }

    #[test]
    fn simple_test_2() {
        let mut tsp: TSP = make_graph_points(&vec![
            (0.0, 0.0).into(),
            (10.0, 0.0).into(),
            (1.0, 0.0).into(),
            (11.0, 0.0).into(),
        ])
        .into();

        assert_eq!(tsp.cost, 29.0);

        tsp.chained_lin_kernighan(1_000, 4);

        assert_eq!(tsp.cost, 11.0);
    }

    #[test]
    fn berlin_chained_lin_kernighan() {
        let mut tsp: TSP = make_graph_points(&berlin()).into();

        tsp.chained_lin_kernighan(1_000, 52);
        // tsp.lin_kernighan(52);
        tsp.update_cost();

        eprintln!("Cost: {}", tsp.cost);
        eprintln!("Cost: {}", tsp.cycle_cost());
        assert!(tsp.cost < 1_000.0);
    }

    #[test]
    fn berlin_held_karp() {
        let mut tsp: TSP = make_graph_points(&berlin()).into();

        tsp.held_karp();

        eprintln!("Cost: {}", tsp.cost);
        eprintln!("Cost: {}", tsp.cycle_cost());
        assert!(tsp.cost < 1_000.0);
    }
}
