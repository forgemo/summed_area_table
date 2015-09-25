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
nalgebra = "0.2.23" // check cargo.toml for dependency for newer version  
```

#### Import the required content.

```rust
use nalgebra::{DMat};
use summed_area_table::{SummedAreaTableSource, SummedAreaTable};
```

#### Creating a summed area table for a 10x10 matrix filled with ones.

```rust
let src: DMat<usize> = DMat::new_ones(10,10);
let table = src.calculate_full_summed_area_table();
```

#### Getting the sum of a specific area

```rust
assert_eq!(100, table.get_sum((0,0),(9,9)));
assert_eq!(50, table.get_sum((0,0),(9,4)));
assert_eq!(25, table.get_sum((0,0),(4,4)));
```

#### Getting the average of a specific area

```rust
assert_eq!(10, table.get_average((0,0),(9,9)));
assert_eq!(10, table.get_average((0,0),(9,4)));
assert_eq!(10, table.get_average((0,0),(4,4)));
```

## Custom Data Source

You can implement the `SummedAreaTableSource`trait for your own types if you want them to have the `calculate_summed_area_table()` method. All you need, is to implement `get_values` and return a DMat<usize> for your data. 

This is how the library does it for `DMat<usize>` type.

```rust
impl SummedAreaTableSource for DMat<usize>{
	fn get_values(&self) -> &DMat<usize> {
		self
	}
}
```
