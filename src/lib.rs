
extern crate nalgebra;
extern crate num;

use nalgebra::{DMat, Indexable};

/// This trait represents the source for a summed area table.
/// Implement this trait for a type to use it as data source for a summed area table.
pub trait SummedAreaTableSource{

	/// This method will be used to access the source data in form af a matrix.
	fn get_values(&self) -> &DMat<usize>;

	/// Calculates and returns the actual summed area table for a given rect.
	/// The arguments 'from' and 'to' represent the rects top-left (inclusive) and bottom-right (inclusive) point.
	fn calculate_summed_area_table(&self, from: (usize, usize), to: (usize, usize)) -> SummedAreaTable{
		let vals = self.get_values();
		let mut table = DMat::new_zeros(vals.nrows(),vals.ncols());

		let (from_x, from_y) = from;
		let (to_x, to_y) = to;
		unsafe {
			for row in (from_y .. to_y+1) {
				for col in (from_x .. to_x+1) {

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


		SummedAreaTable{table: table}
	}

	/// Calculates and returns the actual summed area table for the whole source matrix.
	fn calculate_full_summed_area_table(&self) -> SummedAreaTable{
		let ncols= self.get_values().ncols();
		let nrows= self.get_values().nrows();
		self.calculate_summed_area_table((0,0),(ncols-1, nrows-1))
	}
}

/// This struct represents the result of a summed area table calculation.
pub struct SummedAreaTable {
	pub table: DMat<usize>,
}

impl SummedAreaTable {

	/// Returns the sum for a given area,
	/// that is described by its upper left and lower right point.
	/// It will panic in debug mode if the given points are not within the tables bounds.
	/// It will panic in debug mode if `from` is right of or below `to`.
	/// `from` is a x/y coordinate tuple for the upper left point (inclusive)
	/// `to` is a x/y coordinate tuple for the lower right point (inclusive)
	pub fn get_sum(&self, from: (usize,usize), to: (usize,usize)) -> usize{
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
				sum = sum + self.table.unsafe_at((row1-1,col1-1));
			}
			if col1 > 0 {
				let temp = self.table.unsafe_at((row2,col1-1));

				debug_assert!(temp<=sum, "Overlow-Alarm 1: p1({}/{}) p2({}/{}) temp({}) sum({})",
				col1, row1, col2, row2, temp, sum);

				sum = sum - temp;
			}
			if row1 > 0 {
				let temp = self.table.unsafe_at((row1-1,col2));

				debug_assert!(temp<=sum, "Overlow-Alarm 2: p1({}/{}) p2({}/{}) temp({}) sum({}) ",
				col1, row1, col2, row2, temp, sum);

				sum = sum - temp;
			}
			sum
		}
	}

	/// Returns the average for a given area,
	/// that is described by its upper left and lower right point.
	/// It will panic in debug mode if the given points are not within the tables bounds.
	/// It will panic in debug mode if `from` is right of or below `to`.
	/// `from` is a x/y coordinate tuple for the upper left point (inclusive)
	/// `to` is a x/y coordinate tuple for the lower right point (inclusive)
	pub fn get_average(&self, from: (usize,usize), to: (usize,usize))-> f32{
		let sum = self.get_sum(from,to);
		let data_count = self.get_data_count(from, to);
		sum as f32 / data_count as f32
	}

	/// Returns the number of data points at the given area.
	pub fn get_data_count(&self, from: (usize,usize), to: (usize,usize))-> usize{
		let (from_x, from_y) = from;
		let (to_x, to_y) = to;
		(to_x-from_x+1)*(to_y-from_y+1)
	}

	/// Returns the average for the whole area.
	/// that is described by its upper left and lower right point.
	/// It will panic in debug mode if the given points are not within the tables bounds.
	/// It will panic in debug mode if `from` is right of or below `to`.
	/// `from` is a x/y coordinate tuple for the upper left point (inclusive)
	/// `to` is a x/y coordinate tuple for the lower right point (inclusive)
	pub fn get_overall_average(&self) -> f32{
		self.get_average((0,0),(self.table.ncols()-1,self.table.nrows()-1))
	}

	/// Returns the sum for the whole area.
	/// that is described by its upper left and lower right point.
	/// It will panic in debug mode if the given points are not within the tables bounds.
	/// It will panic in debug mode if `from` is right of or below `to`.
	/// `from` is a x/y coordinate tuple for the upper left point (inclusive)
	/// `to` is a x/y coordinate tuple for the lower right point (inclusive)
	pub fn get_overall_sum(&self) -> usize{
		self.get_sum((0,0),(self.table.ncols()-1,self.table.nrows()-1))
	}

	/// Returns the number of data points at the whole area.
	pub fn get_overall_data_count(&self) -> usize{
		self.table.ncols()*self.table.nrows()
	}
}



impl SummedAreaTableSource for DMat<usize>{
	fn get_values(&self) -> &DMat<usize> {
		self
	}
}



pub mod util {
	use nalgebra::{DMat};
	pub fn vec_to_dmat(vec: &Vec<usize>) -> DMat<usize> {
		DMat::from_col_vec(vec.len(), 1, &vec[..])
	}
}


#[test]
fn zeros() {
	let src: DMat<usize> = DMat::new_zeros(100,100);
	let table = src.calculate_full_summed_area_table();
	assert_eq!(0, table.get_sum((0,0),(99,99)));
}

#[test]
fn ones() {
	let src: DMat<usize> = DMat::from_elem(100,100,1);
	let table = src.calculate_full_summed_area_table();
	assert_eq!(10000, table.get_sum((0,0),(99,99)));
}

#[test]
fn ones_without_first_col_row() {
	let src: DMat<usize> = DMat::from_elem(100,100,1);
	let table = src.calculate_full_summed_area_table();
	assert_eq!(10000-199, table.get_sum((1,1),(99,99)));
}


#[test]
fn twos() {
	let src: DMat<usize> = DMat::from_elem(100,100,2);
	let table = src.calculate_full_summed_area_table();
	assert_eq!(20000, table.get_sum((0,0),(99,99)));
}

#[test]
fn twos_average() {
	let src: DMat<usize> = DMat::from_elem(3,3,2);
	let table = src.calculate_full_summed_area_table();
	assert_eq!(2.0, table.get_average((0,0),(2,2)));
}

#[test]
fn data_count() {
	let src: DMat<usize> = DMat::from_elem(123,321,2);
	let table = src.calculate_full_summed_area_table();
	assert_eq!(123*321, table.get_data_count((0,0),(122,320)));
}

#[test]
fn overall_data_count() {
	let src: DMat<usize> = DMat::from_elem(123,321,2);
	let table = src.calculate_full_summed_area_table();
	assert_eq!(123*321, table.get_overall_data_count());
}

#[test]
fn overall_sum() {
	let src: DMat<usize> = DMat::from_elem(100,100,2);
	let table = src.calculate_full_summed_area_table();
	assert_eq!(20000, table.get_overall_sum());
}

#[test]
fn overall_average() {
	let src: DMat<usize> = DMat::from_elem(3,3,2);
	let table = src.calculate_full_summed_area_table();
	assert_eq!(2.0, table.get_overall_average());
}


#[test]
fn ones_quartered() {
	let src: DMat<usize> = DMat::from_elem(100,100,1);
	let table = src.calculate_full_summed_area_table();
	assert_eq!(2500, table.get_sum((0,0),(49,49)));
	assert_eq!(2500, table.get_sum((50,50),(99,99)));
	assert_eq!(2500, table.get_sum((50,0),(99,49)));
	assert_eq!(2500, table.get_sum((0,50),(49,99)));
}

#[test]
fn twos_quartered() {
	let src: DMat<usize> = DMat::from_elem(100,100,2);
	let table = src.calculate_full_summed_area_table();
	assert_eq!(5000, table.get_sum((0,0),(49,49)));
	assert_eq!(5000, table.get_sum((50,50),(99,99)));
	assert_eq!(5000, table.get_sum((50,0),(99,49)));
	assert_eq!(5000, table.get_sum((0,50),(49,99)));
}

#[test]
fn first_row() {
	let src: DMat<usize> = DMat::from_elem(10,20,1);
	let table = src.calculate_full_summed_area_table();
	assert_eq!(20, table.get_sum((0,0),(19,0)));
}

#[test]
fn first_col() {
	let src: DMat<usize> = DMat::from_elem(50,100,1);
	let table = src.calculate_full_summed_area_table();
	assert_eq!(50, table.get_sum((0,0),(0,49)));
}

#[test]
fn from_to_equal() {
	let src: DMat<usize> = DMat::from_elem(100,100,1);
	let table = src.calculate_full_summed_area_table();
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

	let table = src.calculate_full_summed_area_table();

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
fn vec_to_dmat() {
	let src = util::vec_to_dmat(&vec![0,1,2,3,4,5]);
	assert_eq!(1, src.ncols());
	assert_eq!(6, src.nrows());
	let table = src.calculate_full_summed_area_table();
	assert_eq!(6, table.get_overall_data_count());
	assert_eq!(15, table.get_overall_sum());
}

#[test]
fn src_and_sat_same_size() {
	let src: DMat<usize> = DMat::new_zeros(100,100);
	let table = src.calculate_full_summed_area_table();
	assert_eq!(src.nrows(), table.table.nrows());
	assert_eq!(src.ncols(), table.table.ncols());
}

#[test]
#[should_panic]
fn bound_check_x() {
	let src: DMat<usize> = DMat::new_zeros(50,100);
	let table = src.calculate_full_summed_area_table();
	table.get_sum((0,0),(50,99));
}

#[test]
#[should_panic]
fn bound_check_y() {
	let src: DMat<usize> = DMat::new_zeros(50,100);
	let table = src.calculate_full_summed_area_table();
	table.get_sum((0,0),(49,100));
}

#[test]
#[should_panic]
fn point_order_check1() {
	let src: DMat<usize> = DMat::new_zeros(50,100);
	let table = src.calculate_full_summed_area_table();
	assert_eq!(0, table.get_sum((49,99),(48,98)));
}

#[test]
#[should_panic]
fn point_order_check2() {
	let src: DMat<usize> = DMat::new_zeros(50,100);
	let table = src.calculate_full_summed_area_table();
	assert_eq!(0, table.get_sum((49,99),(48,99)));
}

#[test]
#[should_panic]
fn point_order_check3() {
	let src: DMat<usize> = DMat::new_zeros(50,100);
	let table = src.calculate_full_summed_area_table();
	assert_eq!(0, table.get_sum((49,99),(49,98)));
}
