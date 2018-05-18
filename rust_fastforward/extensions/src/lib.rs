#[macro_use]
extern crate cpython;
extern crate numpy;
extern crate ndarray;
extern crate rand;

use std::f64::INFINITY;

use numpy::*;
use ndarray::prelude::*;
use cpython::{PyResult, Python, PyObject};
use rand::{XorShiftRng};
use rand::distributions::exponential::Exp;
use rand::distributions::IndependentSample;


fn ffs(val: f64, biggerval: f64, otherval: f64, lifetime: f64, _factor: f64, rng: &mut XorShiftRng) -> f64 {
    if val < biggerval && val < otherval { //MIGHT BE WRONG - what if val > biggerval and < otherval
        biggerval + Exp::new(1.0 / lifetime).ind_sample(rng)
    }
    else {
        INFINITY
    }
}

fn fastforward_py(py: Python, transtime: PyArray, nextem: f64, temtime: PyArray, k_trans: f64, factor: f64, _f: PyObject) -> PyResult<PyObject> {
    let mut rng = rand::weak_rng();
    let mut transtime = transtime.as_array_mut().unwrap();
    let temtime = temtime.as_array().unwrap();

    for i in 0.. transtime.shape()[0] {
        for j in 0.. 2usize {
            transtime[[i, j]] = ffs(
                transtime[[i, j]],
                nextem,
                temtime[[i, j]],
                k_trans,
                factor,
                &mut rng
            );
        }
    }

    Ok(py.None())
}

fn deltransgttem_helper_py(py: Python, transtime: PyArray, temtime: PyArray) -> PyResult<usize> {
    let mut transtime: ArrayViewMutD<f64> = transtime.as_array_mut().unwrap();
    let mut temtime: ArrayViewMutD<f64> = temtime.as_array_mut().unwrap();

    let mut len = transtime.shape()[0];
    let mut i = 0;
    while i < len {
        if transtime[[i, 0]] > temtime[[i, 0]] && transtime[[i, 1]] > temtime[[i, 1]] {
            // Delete row, copy row at end into it's place
            transtime[[i, 0]] = transtime[[len - 1, 0]];
            transtime[[i, 1]] = transtime[[len - 1, 1]];
            temtime[[i, 0]] = temtime[[len - 1, 0]];
            temtime[[i, 1]] = temtime[[len - 1, 1]];

            len -= 1;
        }
        else {
            i += 1;
        }
    }

    Ok(len)
}

/* Define module "_rust_ext" */
py_module_initializer!(_rust_fastforward, init_rust_fastforward, PyInit__rust_fastforward, |py, m| {
    m.add(py, "__doc__", "Optimized fast forward code")?;
    m.add(py, "fastforward", py_fn!(py,
        fastforward_py (
            transtime: PyArray,
            nextem: f64,
            temtime: PyArray,
            k_trans: f64,
            factor: f64,
            _f: PyObject
        )))?;
    m.add(py, "deltransgttem_helper", py_fn!(py,
        deltransgttem_helper_py(transtime: PyArray, temtime: PyArray)
    ))?;
    Ok(())
});
