use plotters::prelude::*;
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

fn span_gen(span: (f64, f64), h: f64) -> Vec<f64> {
    let mut l_span = vec![];

    let mut cnt: f64 = span.0;
    while cnt <= span.1 {
        l_span.push(cnt);
        cnt += h;
    }
    l_span.push(cnt);

    l_span
}

// Just analytic
fn analytic<F>(f: &F, x_n: f64, y_n: f64, _h: f64) -> f64
where F:
    Fn(f64, f64) -> f64
{
    f(x_n, y_n)
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
where
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

    // Euler's values
    let ys_e = NumericalSolver::solve(&derivative, &xs, 1., &euler, step);

    // Analytical values
    let ys_a = AnalyticalSolver::solve(&by_hand, &xs, 1., &(), step);

    // Exact Analytical
    let ys_aex = AnalyticalSolver::solve(&by_hand, &xs_exact, 1., &(), 0.05);

    let name = format!("plots/{}.png", number.to_string().replace(".", "_"));



    let root = BitMapBackend::new(&name, (1024, 768)).into_drawing_area();
    root.fill(&WHITE)?;
    let mut chart = ChartBuilder::on(&root)
        .caption("Euler's method and Analytical", ("sans-serif", 50).into_font())
        .margin(5)
        .x_label_area_size(35)
        .y_label_area_size(40)
        .build_cartesian_2d(A..B, 1f64..1.5)?;

    chart.configure_mesh().draw()?;

    // Euler's graph
    chart
        .draw_series(LineSeries::new(
            xs.iter().zip(ys_e.iter()).map(|(&x, &y)| (x, y)),
            &RED
        ))?
        .label("Euler's aproximation")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &RED));

    // Analytical graph
    chart
        .draw_series(LineSeries::new(
            xs_exact.iter().zip(ys_aex.iter()).map(|(&x, &y)| (x, y)),
            &BLUE
        ))?
        .label("Analytical solution")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &BLUE));


    chart
        .configure_series_labels()
        .background_style(&WHITE.mix(0.8))
        .border_style(&BLACK)
        .draw()?;

      root.present()?;

    let e = ys_e.last().unwrap();
    let a = ys_a.last().unwrap();
    println!("Analytical: {}", a);
    println!("Euler: {}", e);
    // Problem with max_error!
    let mx_error = ys_e.iter().zip(ys_a.iter())
                              .map(|(&x, &y)| f64::abs(x - y))
                              .max_by(|a, b| a.partial_cmp(b).unwrap() )
                              .expect("Something went wrong with finding maximum error");

    println!("Max error: {}", mx_error);

    // println!("{:?}", (ys_e, &xs));
    // println!("{:?}", (ys_a, &xs));
    // println!("{:?}", (ys_aex, &xs_exact));

    Ok(())
}
