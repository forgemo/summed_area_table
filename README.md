# summed_area_table
Generic implementation of summed area tables with Rust

You can find more information about summed area tables on Wikipedia.

http://en.wikipedia.org/wiki/Summed_area_table

## Basic usage 

#### Creating a summed area table for a 10x10 matrix filled with ones.

```rust
let src: DMat<usize> = DMat::new_ones(10,10);
let table = src.calculate_summed_area_table();
```

#### Getting the sum of a specific area
```rust
assert_eq!(100, table.get_sum((0,0),(9,9)));
assert_eq!(50, table.get_sum((0,0),(9,4)));
assert_eq!(25, table.get_sum((0,0),(4,4)));
```

## Generic usage

You can implement the `SummedAreaTableSource`trait for your own types if you want them to have the `calculate_summed_area_table()` method. You also need to implement the `SourceValue` trait for the value type of your source. 

This is how the library does it for `DMat<usize>` type. 

```rust
impl SourceValue<usize> for usize {
}

impl SummedAreaTableSource<usize> for DMat<usize>{
	fn get_values(&self) -> &DMat<usize> {
		self
	}
}
```








