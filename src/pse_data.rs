//! atomic weights as found [here ](https://www.ciaaw.org/atomic-weights.htm)
//! If an interval is given in the data source the first value was taken.
//! TODO Check if it is more sound to average over the interval.
pub const PSE_MASSES: [f64; 19] = [
    0.0,               // place holder
    1.007_825_031_898, // h
    4.002_603_241,     // he
    6.997,             // li
    9.012_183_1,       // be
    10.806,            // b
    12.009_6,          // c
    14.006_43,         // n
    15.999_03,         // o
    18.998_403_162,    // f
    20.1797,           // ne
    22.989_769_282,    // na
    24.304,            // mg
    26.981_538_4,      // al
    28.084,            // si
    30.973_761_998,    // p
    32.059,            // s
    35.446,            // cl
    39.792,            // ar
];
