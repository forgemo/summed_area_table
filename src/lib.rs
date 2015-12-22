extern crate nalgebra;

use nalgebra::{DMat};

pub trait SourceValue : Copy {
	fn as_f64(self) -> f64;
}

impl SourceValue for u8 { fn as_f64(self) -> f64 { self as f64 } }
impl SourceValue for i8 { fn as_f64(self) -> f64 { self as f64 } }
impl SourceValue for u32 { fn as_f64(self) -> f64 { self as f64 } }
impl SourceValue for i32 { fn as_f64(self) -> f64 { self as f64 } }
impl SourceValue for u64 { fn as_f64(self) -> f64 { self as f64 } }
impl SourceValue for i64 { fn as_f64(self) -> f64 { self as f64 } }
impl SourceValue for usize { fn as_f64(self) -> f64 { self as f64 } }
impl SourceValue for isize { fn as_f64(self) -> f64 { self as f64 } }
impl SourceValue for f32 { fn as_f64(self) -> f64 { self as f64 } }
impl SourceValue for f64 { fn as_f64(self) -> f64 { self } }

impl <'a>SourceValue for &'a u8 { fn as_f64(self) -> f64 { *self as f64 } }
impl <'a>SourceValue for &'a i8 { fn as_f64(self) -> f64 { *self as f64 } }
impl <'a>SourceValue for &'a u32 { fn as_f64(self) -> f64 { *self as f64 } }
impl <'a>SourceValue for &'a i32 { fn as_f64(self) -> f64 { *self as f64 } }
impl <'a>SourceValue for &'a u64 { fn as_f64(self) -> f64 { *self as f64 } }
impl <'a>SourceValue for &'a i64 { fn as_f64(self) -> f64 { *self as f64 } }
impl <'a>SourceValue for &'a usize { fn as_f64(self) -> f64 { *self as f64 } }
impl <'a>SourceValue for &'a isize { fn as_f64(self) -> f64 { *self as f64 } }
impl <'a>SourceValue for &'a f32 { fn as_f64(self) -> f64 { *self as f64 } }
impl <'a>SourceValue for &'a f64 { fn as_f64(self) -> f64 { *self } }

/// This trait represents the source for a summed area table.
/// Implement this trait for a type to use it as data source for a summed area table.
pub trait SummedAreaTableSource<T: SourceValue>{

	/// This method will be used to access the source data in form af a matrix.
	fn at(&self, x:usize, y:usize) -> &T;
	fn height(&self) -> usize;
	fn width(&self) -> usize;

	/// Calculates and returns the actual summed area table for a given rect.
	/// The arguments 'from' and 'to' represent the rects top-left (inclusive) and bottom-right (inclusive) point.
	fn calculate_summed_area_table(&self, from: (usize, usize), to: (usize, usize)) -> SummedAreaTable{
		let mut table:DMat<f64> = DMat::new_zeros(self.height(),self.width());

		let (from_x, from_y) = from;
		let (to_x, to_y) = to;

		for row in from_y .. to_y+1 {
			for col in from_x .. to_x+1 {

				let mut sum = self.at(col, row).as_f64();

				if row>0 {
					sum = sum + table[(row-1, col)];
				}
				if col>0 {
					sum = sum+ table[(row, col-1)];
				}
				if row>0 && col>0 {
					sum = sum - table[(row-1, col-1)];
				}
				table[(row,col)] = sum;
			}
		}

		SummedAreaTable{table: table}
	}

	/// Calculates and returns the actual summed area table for the whole source matrix.
	fn calculate_full_summed_area_table(&self) -> SummedAreaTable{
		let ncols= self.width();
		let nrows= self.height();
		self.calculate_summed_area_table((0,0),(ncols-1, nrows-1))
	}
}

/// This struct represents the result of a summed area table calculation.
pub struct SummedAreaTable {
	pub table: DMat<f64>,
}

impl SummedAreaTable {

	/// Returns the sum for a given area,
	/// that is described by its upper left and lower right point.
	/// It will panic in debug mode if the given points are not within the tables bounds.
	/// It will panic in debug mode if `from` is right of or below `to`.
	/// `from` is a x/y coordinate tuple for the upper left point (inclusive)
	/// `to` is a x/y coordinate tuple for the lower right point (inclusive)
	pub fn get_sum(&self, from: (usize,usize), to: (usize,usize)) -> f64{
		let (col1, row1) = from;
		let (col2, row2) = to;

		debug_assert!(row1 <= row2 && col1 <= col2, "`from` ({}/{}) must not be right of or below `to`({}/{})", col1, row1, col2, row2);

		debug_assert!( {
			let ncols = self.table.ncols();
			let nrows = self.table.nrows();
			col1 < ncols && col2 < ncols && row1 < nrows && row2 < nrows
		},"`from` ({}/{}) or `to` ({}/{}) not within table bounds [(0/0)..({}/{})]", col1, row1, col2, row2,  self.table.ncols()-1,  self.table.nrows()-1);


		let mut sum = self.table[(row2,col2)];

		if col1 > 0 && row1 > 0 {
			sum = sum + self.table[(row1-1,col1-1)];
		}
		if col1 > 0 {
			let temp = self.table[(row2,col1-1)];

			debug_assert!(temp<=sum, "Overlow-Alarm 1: p1({}/{}) p2({}/{}) temp({}) sum({})",
			col1, row1, col2, row2, temp, sum);

			sum = sum - temp;
		}
		if row1 > 0 {
			let temp = self.table[(row1-1,col2)];

			debug_assert!(temp<=sum, "Overlow-Alarm 2: p1({}/{}) p2({}/{}) temp({}) sum({}) ",
			col1, row1, col2, row2, temp, sum);

			sum = sum - temp;
		}
		sum
	}


	/// Returns the average for a given area,
	/// that is described by its upper left and lower right point.
	/// It will panic in debug mode if the given points are not within the tables bounds.
	/// It will panic in debug mode if `from` is right of or below `to`.
	/// `from` is a x/y coordinate tuple for the upper left point (inclusive)
	/// `to` is a x/y coordinate tuple for the lower right point (inclusive)
	pub fn get_average(&self, from: (usize,usize), to: (usize,usize))-> f64{
		let sum = self.get_sum(from,to);
		let data_count = self.get_data_count(from, to);
		sum/data_count as f64
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
	pub fn get_overall_average(&self) -> f64{
		self.get_average((0,0),(self.table.ncols()-1,self.table.nrows()-1))
	}

	/// Returns the sum for the whole area.
	/// that is described by its upper left and lower right point.
	/// It will panic in debug mode if the given points are not within the tables bounds.
	/// It will panic in debug mode if `from` is right of or below `to`.
	/// `from` is a x/y coordinate tuple for the upper left point (inclusive)
	/// `to` is a x/y coordinate tuple for the lower right point (inclusive)
	pub fn get_overall_sum(&self) -> f64{
		self.get_sum((0,0),(self.table.ncols()-1,self.table.nrows()-1))
	}

	/// Returns the number of data points at the whole area.
	pub fn get_overall_data_count(&self) -> usize{
		self.table.ncols()*self.table.nrows()
	}
}



impl <T: SourceValue>SummedAreaTableSource<T> for DMat<T>{
	fn at(&self, x: usize, y: usize) -> &T {
		&self[(y,x)]
	}
	fn height(&self) -> usize {
		self.nrows()
	}
	fn width(&self) -> usize {
		self.ncols()
	}
}
pub struct VecSource<'a,T:'a>{
	vec: &'a [T],
	height: usize,
	width: usize
}
impl <'a, T: SourceValue>VecSource<'a, T>{
	pub fn new(vec: &'a [T], width: usize, height: usize) -> VecSource<'a, T> {
		VecSource{vec: vec, width: width, height: height}
	}
}
impl <'a, T: SourceValue>SummedAreaTableSource<T> for VecSource<'a, T>{

	fn at(&self, x: usize, y: usize) -> &T {
		&self.vec[util::map_2d_to_1d(x, y, self.width)]
	}
	fn height(&self) -> usize {
		self.height
	}
	fn width(&self) -> usize {
		self.width
	}
}

pub mod util {
	use nalgebra::{DMat};
	pub fn vec_to_dmat(vec: &Vec<usize>) -> DMat<usize> {
		DMat::from_col_vec(vec.len(), 1, &vec[..])
	}
	pub fn map_2d_to_1d(x: usize, y: usize, width: usize)->usize{
	    y * width + x
	}
}


#[test]
fn zeros() {
	let src: DMat<usize> = DMat::new_zeros(100,100);
	let table = src.calculate_full_summed_area_table();
	assert_eq!(0.0, table.get_sum((0,0),(99,99)));
}

#[test]
fn ones() {
	let src: DMat<usize> = DMat::from_elem(100,100,1);
	let table = src.calculate_full_summed_area_table();
	assert_eq!(10000.0, table.get_sum((0,0),(99,99)));
}

#[test]
fn ones_without_first_col_row() {
	let src: DMat<usize> = DMat::from_elem(100,100,1);
	let table = src.calculate_full_summed_area_table();
	assert_eq!(10000.0-199.0, table.get_sum((1,1),(99,99)));
}


#[test]
fn twos() {
	let src: DMat<usize> = DMat::from_elem(100,100,2);
	let table = src.calculate_full_summed_area_table();
	assert_eq!(20000.0, table.get_sum((0,0),(99,99)));
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
	assert_eq!(20000.0, table.get_overall_sum());
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
	assert_eq!(2500.0, table.get_sum((0,0),(49,49)));
	assert_eq!(2500.0, table.get_sum((50,50),(99,99)));
	assert_eq!(2500.0, table.get_sum((50,0),(99,49)));
	assert_eq!(2500.0, table.get_sum((0,50),(49,99)));
}

#[test]
fn ones_quartered_vec() {
	let mat = DMat::from_elem(100,100,1);
	let src = VecSource::new(mat.as_vec(), 100,100);
	let table = src.calculate_full_summed_area_table();
	assert_eq!(2500.0, table.get_sum((0,0),(49,49)));
	assert_eq!(2500.0, table.get_sum((50,50),(99,99)));
	assert_eq!(2500.0, table.get_sum((50,0),(99,49)));
	assert_eq!(2500.0, table.get_sum((0,50),(49,99)));
}

#[test]
fn twos_quartered() {
	let src: DMat<usize> = DMat::from_elem(100,100,2);
	let table = src.calculate_full_summed_area_table();
	assert_eq!(5000.0, table.get_sum((0,0),(49,49)));
	assert_eq!(5000.0, table.get_sum((50,50),(99,99)));
	assert_eq!(5000.0, table.get_sum((50,0),(99,49)));
	assert_eq!(5000.0, table.get_sum((0,50),(49,99)));
}

#[test]
fn first_row() {
	let src: DMat<usize> = DMat::from_elem(10,20,1);
	let table = src.calculate_full_summed_area_table();
	assert_eq!(20.0, table.get_sum((0,0),(19,0)));
}

#[test]
fn first_col() {
	let src: DMat<usize> = DMat::from_elem(50,100,1);
	let table = src.calculate_full_summed_area_table();
	assert_eq!(50.0, table.get_sum((0,0),(0,49)));
}

#[test]
fn from_to_equal() {
	let src: DMat<usize> = DMat::from_elem(100,100,1);
	let table = src.calculate_full_summed_area_table();
	assert_eq!(1.0, table.get_sum((0,0),(0,0)));
	assert_eq!(1.0, table.get_sum((50,50),(50,50)));
	assert_eq!(1.0, table.get_sum((99,99),(99,99)));
}

#[test]
fn custom() {

	let src: DMat<f64> = DMat::from_row_vec(5,5, &[
		5.0,2.0,3.0,4.0,1.0,
		1.0,5.0,4.0,2.0,3.0,
		2.0,2.0,1.0,3.0,4.0,
		3.0,5.0,6.0,4.0,5.0,
		4.0,1.0,3.0,2.0,6.0
	]);
	let expected_table: DMat<f64> = DMat::from_row_vec(5,5, &[
		5.0,7.0,10.0,14.0,15.0,
		6.0,13.0,20.0,26.0,30.0,
		8.0,17.0,25.0,34.0,42.0,
		11.0,25.0,39.0,52.0,65.0,
		15.0,30.0,47.0,62.0,81.0
	]);

	let table = src.calculate_full_summed_area_table();

	for row in 0 .. 5 {
		for col in 0 .. 5 {
			assert_eq!(expected_table[(row,col)], table.table[(row,col)]);
		}
	};
	for x in 0 .. 5 {
		for y in 0 .. 5 {
			assert_eq!(expected_table[(y,x)].as_f64(), table.get_sum((0,0), (x,y)));
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
	assert_eq!(15.0, table.get_overall_sum());
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
	table.get_sum((49,99),(48,98));
}

#[test]
#[should_panic]
fn point_order_check2() {
	let src: DMat<usize> = DMat::new_zeros(50,100);
	let table = src.calculate_full_summed_area_table();
	table.get_sum((49,99),(48,99));
}

#[test]
#[should_panic]
fn point_order_check3() {
	let src: DMat<usize> = DMat::new_zeros(50,100);
	let table = src.calculate_full_summed_area_table();
	table.get_sum((49,99),(49,98));
}
