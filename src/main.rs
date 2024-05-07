use plotters::prelude::*;
// 4xyy' - 3y^2 + x^2 = 0
// =>
// y' = (3y/4) - (x^2)/4(y)
//
// y(1) = 1
// x \in [1;2]

// rhs y' = RHS
// for numerical
fn derivative(x: f64, y: f64) -> f64 {
    ((3. * y) / (4. * x)) - (x / (4. * y))
}

// y(x) = sqrt(-x^2 + c_1*x^(3/2))
// where c_1 = 2
fn by_hand(x: f64, y: f64) -> f64 {
    (-(x.powf(2.)) + 2.*x.powf(3./2.)).sqrt()
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


// consider is something like 0.01 (1.0 -> 1.01 -> 1.02 -> ...)
// WHY X_SPAN is (i32, i32)?
// It's all beacuse of float vals arithmetic which not allowing provide
// nice api for this
fn solve_ivp<F1, F2>(f: &F1, x_span: (f64, f64), y_0: f64, h: f64, solver: &F2) -> (Vec<f64>, Vec<f64>)
    where
    F1: Fn(f64, f64) -> f64,
    F2: Fn(&F1, f64, f64, f64) -> f64,
{
    // let xs = ((x_span.0)..=(x_span.1)).map(|x| x as f64 * h);
    // let xs_cpy: Vec<f64> = xs.clone().collect();
    let mut xs = vec![];

    let mut cnt: f64 = x_span.0;
    while cnt < x_span.1 {
        xs.push(cnt);
        cnt += h;
    }
    xs.push(cnt);

    let mut ys = vec![0.; xs.len()];
    ys[0] = y_0;

    for i in 0..xs.len() - 1 {
        ys[i + 1] = solver(f, xs[i], ys[i], h);
    }


    (ys, xs)
}


fn main() -> Result<(), Box<dyn std::error::Error>>{
    const STEP: f64 = 0.25;
    let name = format!("plots/{}.png", STEP.to_string().replace(".", "_"));
    // Euler's values
    let (ys_e, xs_e) = solve_ivp(&derivative, (1.0, 2.0), 1., STEP, &euler);
    // Analytical values
    let (ys_a, xs_a) = solve_ivp(&by_hand, (1.0, 2.0), 1., 0.001, &analytic);

    let root = BitMapBackend::new(&name, (1024, 768)).into_drawing_area();
    root.fill(&WHITE)?;
    let mut chart = ChartBuilder::on(&root)
        .caption("Euler's method and Analytical", ("sans-serif", 50).into_font())
        .margin(5)
        .x_label_area_size(35)
        .y_label_area_size(40)
        .build_cartesian_2d(1f64..2f64, 0f64..3f64)?;

    chart.configure_mesh().draw()?;

    // Euler's graph
    chart
        .draw_series(LineSeries::new(
            xs_e.iter().zip(ys_e.iter()).map(|(&x, &y)| (x, y)),
            &RED
        ))?
        .label("Euler's aproximation")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &RED));

    // Analytical graph
    chart
        .draw_series(LineSeries::new(
            xs_a.iter().zip(ys_a.iter()).map(|(&x, &y)| (x, y)),
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

    println!("Analytical: {}", ys_a.last().unwrap());
    println!("Euler: {}", ys_e.last().unwrap());

    Ok(())
}
