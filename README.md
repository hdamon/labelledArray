# labelledArray

A Matlab object class for manipulating structured arrays.

The primary goal of the labelledArray class is to provide data labelling and access capabilities similar to the pandas or xarray packages available for Python.

labelledArray equips numeric arrays with four additional properties:
1) dimNames  : Each dimension of the array can be named
2) dimLabels : Each element of each dimension can be individually labelled
3) dimUnits  : Each dimension can have a physical unit associated with it
4) dimValues : Each dimension can have a secondary set of values associated with it.

labelledArrays can be referenced using relatively arbitrary combinations of dimension names, labels, and indices, to allow more natural syntax when applying operations to relatively complex structures.

labelledArray overloads the most common mathematical operators, and implementing implicit dimension expansion through the use of bsxfun().

Currently implemented functions include:
    permute()
    cat()
    plus()
    minus()
    ldivide()
    rdivide()
    times()
    mean()
    std()

Helper methods are included to easily add additional functionality.



