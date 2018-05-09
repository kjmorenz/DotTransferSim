#[macro_use]
extern crate cpython;
extern crate numpy;
extern crate ndarray;
extern crate rand;

use std::f64::INFINITY;

// use numpy::*;
use cpython::{PyResult, Python, PyObject, PyList};
use rand::XorShiftRng;
use rand::distributions::exponential::Exp;
use rand::distributions::IndependentSample;
use cpython::ObjectProtocol;
use cpython::ToPyObject;
use cpython::PythonObject;


fn ffs(val: f64, biggerval: f64, otherval: f64, lifetime: f64, _factor: f64, rng: &mut XorShiftRng) -> f64 {
    if val < biggerval && val < otherval {
        biggerval + Exp::new(1.0 / lifetime).ind_sample(rng)
    }
    else {
        INFINITY
    }
}

fn fastforward_py(py: Python, transtime: PyList, nextem: f64, temtime: PyList, k_trans: f64, factor: f64, _f: PyObject) -> PyResult<PyObject> {
    let mut rng = rand::weak_rng();

    for i in 0.. transtime.len(py) {
        for j in 0.. 2usize {
            let transtime_i: PyList = transtime.get_item(py, i).extract(py)?;
            transtime_i.set_item(py, j, ffs(
                transtime_i.get_item(py, j).extract(py)?,
                nextem,
                temtime.get_item(py, i).get_item(py, j)?.extract(py)?,
                k_trans,
                factor,
                &mut rng
            ).to_py_object(py).into_object());
        }
    }

    Ok(py.None())
}


/* Define module "_rust_ext" */
py_module_initializer!(_rust_fastforward, init_rust_fastforward, PyInit__rust_fastforward, |py, m| {
    m.add(py, "__doc__", "Optimized fast forward code")?;
    m.add(py, "fastforward", py_fn!(py,
        fastforward_py (
            transtime: PyList,
            nextem: f64,
            temtime: PyList,
            k_trans: f64,
            factor: f64,
            _f: PyObject
        )))?;
    Ok(())
});
