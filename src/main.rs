use plotly::{
    common::{Mode, Line, DashType},
    layout::{Axis, Layout},
    Plot, Scatter
};
use std::io;
const A: f64 = 1.0;
const B: f64 = 2.0;


// 4xyy' - 3y^2 + x^2 = 0
// =>
// y' = ...
// for numerical
fn derivative(x: f64, y: f64) -> f64 {
    ((3. * y) / (4. * x)) - (x / (4. * y))
}

// y(x) = sqrt(-x^2 + c_1*x^(3/2))
// where c_1 = 2
fn by_hand(x: f64) -> f64 {
    (-x.powf(2.) + 2.*x.powf(3./2.)).sqrt()
}

// Generate range
fn span_gen(span: (f64, f64), h: f64) -> Vec<f64> {
  let mut l_span = vec![];
  let mut cnt: f64 = span.0;
  let tolerance = 1e-6; // Adjust tolerance as needed

  while cnt < span.1 + tolerance {
    l_span.push(cnt);
    cnt += h;
  }
  l_span
}

// Euler's method
fn euler<F>(f: &F, x_n: f64, y_n: f64, h: f64) -> f64
where F:
    Fn(f64, f64) -> f64
{
    y_n + h * f(x_n, y_n)
}

/// F1 is a function itself
/// F2 is a solver function which applies F1
trait Solver<F1, F2>
{
    fn solve(f: &F1, xs: &Vec<f64>, y_0: f64, solver: &F2, h: f64) -> Vec<f64>;
}

struct NumericalSolver;

impl<F1, F2> Solver<F1, F2> for NumericalSolver
where F1: Fn(f64, f64) -> f64,
      F2: Fn(&F1, f64, f64, f64) -> f64
{
    fn solve(f: &F1, xs: &Vec<f64>, y_0: f64, solver: &F2, h: f64) -> Vec<f64> {
        let mut ys = vec![0.; xs.len()];
        ys[0] = y_0;

        for i in 0..xs.len() - 1 {
            ys[i + 1] = solver(f, xs[i], ys[i], h);
        }

        ys
    }
}

struct AnalyticalSolver;

impl<F1> Solver<F1, ()> for AnalyticalSolver
where F1: Fn(f64) -> f64
{
    fn solve(f: &F1, xs: &Vec<f64>, _y_0: f64, _solver: &(), _h: f64) -> Vec<f64> {
        let mut ys = vec![0.; xs.len()];

        for i in 0..xs.len() {
            ys[i] = f(xs[i]);
        }

        ys
    }
}


fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("Please enter a step size:");
    let mut input = String::new();
    io::stdin().read_line(&mut input)
        .expect("Failed to read line");

    let number: i32 = input.trim().parse()
        .expect("Please enter a valid number");

    let step: f64 = (B - A) / (number as f64);


    let xs = span_gen((A, B), step);
    let xs_exact = span_gen((A, B), 0.05);
    // println!("xs_e: {:#?}", xs_exact);

    // Euler's values
    let ys_e = NumericalSolver::solve(&derivative, &xs, 1., &euler, step);

    // Analytical values
    let ys_a = AnalyticalSolver::solve(&by_hand, &xs, 1., &(), step);

    // Exact Analytical
    let ys_aex = AnalyticalSolver::solve(&by_hand, &xs_exact, 1., &(), 0.05);

    let mut plot = Plot::new();
    let trace_euler = Scatter::new(xs.clone(), ys_e.clone())
        .mode(Mode::LinesMarkers)
        .line(Line::new().dash(DashType::Dash))
        .name("Euler's method");

    let trace_analytical = Scatter::new(xs_exact.clone(), ys_aex.clone())
        .mode(Mode::Lines)
        .name("Analitycal method");

    let layout = Layout::new()
        .title("Euler's aproximation".into())
        .x_axis(Axis::new().title("x".into()).range(vec![0., 2.1]))
        .y_axis(Axis::new().title("y".into()).range(vec![0.9, 2.]))
        .height(600)
        .width(1200);

    plot.set_layout(layout);
    plot.add_trace(trace_euler);
    plot.add_trace(trace_analytical);

    plot.show();

    let mx_error = ys_e.iter().zip(ys_a.iter())
                              .map(|(&x, &y)| f64::abs(x - y))
                              .max_by(|a, b| a.partial_cmp(b).unwrap() )
                              .expect("Something went wrong with finding maximum error");

    println!("Max error: {}", mx_error);

    Ok(())
}
