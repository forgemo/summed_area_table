#![feature(test)]

extern crate nalgebra;
extern crate test;
extern crate num;

use num::Zero;
use nalgebra::{DMat, Indexable};
use std::ops::{Add, Sub};
use std::fmt::{Display};

/// This trait describes the type of values within a summed area table source.
/// Implementing this trait for your type allows it to become a value within a source.
pub trait SourceValue<T>: Display  + Clone + Copy + Zero + Add<T,Output=T> + Sub<T,Output=T> {
}


/// This trait describes a generic source for a summed area table.
/// Implementing this trait for your type allows it to become the source for a summed area table calculation.
pub trait SummedAreaTableSource<T: SourceValue<T>>{

	/// Returns a matrix of source values of type T.
	/// The calculated summed area table will be based on these values.
	fn get_values(&self) -> &DMat<T>;

	/// Calculates and returns the summed area table for the source matrix.
	fn calculate_summed_area_table(&self) -> SummedAreaTable<T>{
		let vals = self.get_values();
		let mut table = DMat::new_zeros(vals.nrows(),vals.ncols());

		unsafe {
			for row in (0 ..vals.nrows()) {
				for col in (0 ..vals.ncols()) {

					let mut sum = vals.unsafe_at((row, col));

					if row>0 {
						sum = sum+ table.unsafe_at((row-1, col));
					}
					if col>0 {
						sum = sum+ table.unsafe_at((row, col-1));
					}
					if row>0 && col>0 {
						sum = sum - table.unsafe_at((row-1, col-1));
					}
					table.unsafe_set((row,col), sum);
				}
			}
		}


		SummedAreaTable::<T>{table: table}
	}
}

/// This struct is the result of a summed area table calculation.
pub struct SummedAreaTable<T: SourceValue<T>> {
	pub table: DMat<T>,
}

impl <T: SourceValue<T>>SummedAreaTable<T> {

	/// Returns the sum for a given area,
	/// that is described by its upper left and lower right point.
	/// It will panic in debug mode if the given points are not within the tables bounds.
	/// It will panic in debug mode if `from` is right of or below `to`.
	/// `from` is a x/y coordinate tuple for the upper left point (inclusive)
	/// `to` is a x/y coordinate tuple for the lower right point (inclusive)
	pub fn get_sum(&self, from: (usize,usize), to: (usize,usize)) -> T{
		let (col1, row1) = from;
		let (col2, row2) = to;

		debug_assert!(row1 <= row2 && col1 <= col2, "`from` ({}/{}) must not be right of or below `to`({}/{})", col1, row1, col2, row2);

		debug_assert!( {
			let ncols = self.table.ncols();
			let nrows = self.table.nrows();
			col1 < ncols && col2 < ncols && row1 < nrows && row2 < nrows
		},"`from` ({}/{}) or `to` ({}/{}) not within table bounds [(0/0)..({}/{})]", col1, row1, col2, row2,  self.table.ncols()-1,  self.table.nrows()-1);

		unsafe {
			let mut sum = self.table.unsafe_at((row2,col2));

			if col1 > 0 && row1 > 0 {
				sum = sum + self.table.unsafe_at((col1-1,row1-1));
			}
			if col1 > 0 {
				sum = sum - self.table.unsafe_at((col1-1,row2));
			}
			if row1 > 0 {
				sum = sum - self.table.unsafe_at((col2,row1-1));
			}
			sum
		}
	}
}


impl SourceValue<usize> for usize {
}

impl SummedAreaTableSource<usize> for DMat<usize>{
	fn get_values(&self) -> &DMat<usize> {
		self
	}
}


#[test]
fn zeros() {
	let src: DMat<usize> = DMat::new_zeros(100,100);
	let table = src.calculate_summed_area_table();
	assert_eq!(0, table.get_sum((0,0),(99,99)));
}

#[test]
fn ones() {
	let src: DMat<usize> = DMat::from_elem(100,100,1);
	let table = src.calculate_summed_area_table();
	assert_eq!(10000, table.get_sum((0,0),(99,99)));
}

#[test]
fn ones_without_first_col_row() {
	let src: DMat<usize> = DMat::from_elem(100,100,1);
	let table = src.calculate_summed_area_table();
	assert_eq!(10000-199, table.get_sum((1,1),(99,99)));
}


#[test]
fn twos() {
	let src: DMat<usize> = DMat::from_elem(100,100,2);
	let table = src.calculate_summed_area_table();
	assert_eq!(20000, table.get_sum((0,0),(99,99)));
}


#[test]
fn ones_quartered() {
	let src: DMat<usize> = DMat::from_elem(100,100,1);
	let table = src.calculate_summed_area_table();
	assert_eq!(2500, table.get_sum((0,0),(49,49)));
	assert_eq!(2500, table.get_sum((50,50),(99,99)));
	assert_eq!(2500, table.get_sum((50,0),(99,49)));
	assert_eq!(2500, table.get_sum((0,50),(49,99)));
}

#[test]
fn twos_quartered() {
	let src: DMat<usize> = DMat::from_elem(100,100,2);
	let table = src.calculate_summed_area_table();
	assert_eq!(5000, table.get_sum((0,0),(49,49)));
	assert_eq!(5000, table.get_sum((50,50),(99,99)));
	assert_eq!(5000, table.get_sum((50,0),(99,49)));
	assert_eq!(5000, table.get_sum((0,50),(49,99)));
}

#[test]
fn first_row() {
	let src: DMat<usize> = DMat::from_elem(10,20,1);
	let table = src.calculate_summed_area_table();
	assert_eq!(20, table.get_sum((0,0),(19,0)));
}

#[test]
fn first_col() {
	let src: DMat<usize> = DMat::from_elem(50,100,1);
	let table = src.calculate_summed_area_table();
	assert_eq!(50, table.get_sum((0,0),(0,49)));
}

#[test]
fn from_to_equal() {
	let src: DMat<usize> = DMat::from_elem(100,100,1);
	let table = src.calculate_summed_area_table();
	assert_eq!(1, table.get_sum((0,0),(0,0)));
	assert_eq!(1, table.get_sum((50,50),(50,50)));
	assert_eq!(1, table.get_sum((99,99),(99,99)));
}

#[test]
fn custom() {

	let src = DMat::from_row_vec(5,5, &[
		5,2,3,4,1,
		1,5,4,2,3,
		2,2,1,3,4,
		3,5,6,4,5,
		4,1,3,2,6
	]);
	let expected_table = DMat::from_row_vec(5,5, &[
		5,7,10,14,15,
		6,13,20,26,30,
		8,17,25,34,42,
		11,25,39,52,65,
		15,30,47,62,81
	]);

	let table = src.calculate_summed_area_table();

	unsafe{
		for row in 0 .. 5 {
			for col in 0 .. 5 {
				assert_eq!(expected_table.unsafe_at((row,col)), table.table.unsafe_at((row,col)));
			}
		};
		for x in 0 .. 5 {
			for y in 0 .. 5 {
				assert_eq!(expected_table.unsafe_at((y,x)), table.get_sum((0,0), (x,y)));
			}
		}
	}
}

#[test]
#[should_panic]
fn bound_check_x() {
	let src: DMat<usize> = DMat::new_zeros(50,100);
	let table = src.calculate_summed_area_table();
	table.get_sum((0,0),(50,99));
}

#[test]
#[should_panic]
fn bound_check_y() {
	let src: DMat<usize> = DMat::new_zeros(50,100);
	let table = src.calculate_summed_area_table();
	table.get_sum((0,0),(49,100));
}

#[test]
#[should_panic]
fn point_order_check1() {
	let src: DMat<usize> = DMat::new_zeros(50,100);
	let table = src.calculate_summed_area_table();
	assert_eq!(0, table.get_sum((49,99),(48,98)));
}

#[test]
#[should_panic]
fn point_order_check2() {
	let src: DMat<usize> = DMat::new_zeros(50,100);
	let table = src.calculate_summed_area_table();
	assert_eq!(0, table.get_sum((49,99),(48,99)));
}

#[test]
#[should_panic]
fn point_order_check3() {
	let src: DMat<usize> = DMat::new_zeros(50,100);
	let table = src.calculate_summed_area_table();
	assert_eq!(0, table.get_sum((49,99),(49,98)));
}

#[bench]
fn large_1k_matrix(b: &mut test::Bencher) {
	let src: DMat<usize> = DMat::from_elem(1000,1000,1);
	b.iter(|| {
		src.calculate_summed_area_table();
	});
}

#[bench]
fn large_2k_matrix(b: &mut test::Bencher) {
	let src: DMat<usize> = DMat::from_elem(2000,2000,1);
	b.iter(|| {
		src.calculate_summed_area_table();
	})
}

#[bench]
fn large_4k_matrix(b: &mut test::Bencher) {
	let src: DMat<usize> = DMat::from_elem(4000,4000,1);
	b.iter(|| {
		src.calculate_summed_area_table();
	})
}
