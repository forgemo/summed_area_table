extern crate nalgebra;

use nalgebra::{DMat, Zero, Indexable};
use std::ops::{Add, Sub};
use std::fmt::{Display};

pub trait SourceValue<T>: Display + Zero + Clone + Copy + Add<T,Output=T> + Sub<T,Output=T> {
}

pub trait SummedAreaTableSource<T: SourceValue<T>>{

	fn get_values(&self) -> &DMat<T>;

	fn calculate_summed_area_table(&self) -> SummedAreaTable<T>{
		let vals = self.get_values();
		let mut table = DMat::new_zeros(vals.nrows(),vals.ncols());

		for row in (0 ..vals.nrows()) {
			for col in (0 ..vals.ncols()) {
				let mut sum = vals.at((row, col));

				if row>0 {
					sum = sum+ table.at((row-1, col));
				}
				if col>0 {
					sum = sum+ table.at((row, col-1));
				}
				if row>0 && col>0 {
					sum = sum - table.at((row-1, col-1));
				}
				table.set((row,col), sum);
			}
		}

		SummedAreaTable::<T>{table: table}
	}
}

pub struct SummedAreaTable<T: SourceValue<T>> {
	table: DMat<T>,
}

impl <T: SourceValue<T>>SummedAreaTable<T> {
	fn get_sum(&self, from: (usize,usize), to: (usize,usize)) -> T{
		let (x1, y1) = from;
		let (x2, y2) = to;

		let mut sum = self.table.at((x2,y2));

		if x1 > 0 {
			sum = sum - self.table.at((x1-1,y2));
		}
		if y1 > 0 {
			sum = sum - self.table.at((x2,y1-1));
		}
		if x1 > 0 && y1 > 0 {
			sum = sum + self.table.at((x1-1,y1-1));
		}
		sum
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

