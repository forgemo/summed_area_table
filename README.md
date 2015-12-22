[![Build Status](https://travis-ci.org/forgemo/summed_area_table.svg?branch=master)](https://travis-ci.org/forgemo/summed_area_table)

# summed_area_table
Generic implementation of summed area tables with Rust

You can find more information about summed area tables on Wikipedia.

http://en.wikipedia.org/wiki/Summed_area_table

## Basic usage

#### Dependecies for your cargo.toml

```rust
[dependencies]
summed-area-table = "*"
nalgebra = "0.4.0" // this is not needed if you don't want to use nalgebra::DMat
```

#### Import the required content.

```rust
use nalgebra::{DMat}; // this is not needed if you don't want to use nalgebra::DMat
use summed_area_table::{SummedAreaTableSource, SummedAreaTable};
```

#### Creating a summed area table for a 10x10 matrix filled with ones.

```rust
let src: DMat<usize> = DMat::new_ones(10,10);
let table = src.calculate_full_summed_area_table();
```

#### Getting the sum of a specific area

```rust
assert_eq!(100.0, table.get_sum((0,0),(9,9)));
assert_eq!(50.0, table.get_sum((0,0),(9,4)));
assert_eq!(25.0, table.get_sum((0,0),(4,4)));
```

#### Getting the average of a specific area

```rust
assert_eq!(10.0, table.get_average((0,0),(9,9)));
assert_eq!(10.0, table.get_average((0,0),(9,4)));
assert_eq!(10.0, table.get_average((0,0),(4,4)));
```

## Custom Data Source

You can implement the `SummedAreaTableSource` trait for your own types if you want them to have the `calculate_summed_area_table()` method. All you need, is to implement the `SummedAreaTableSource<T>` trait.
`T` need to implement the `SourceValue` trait. This library, however, comes with implementations for all primitive numeric types.

The following code shows how this library implements `SummedAreaTableSource<T>` for the `nalgebra::DMat<T>` type.

```rust
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
```
