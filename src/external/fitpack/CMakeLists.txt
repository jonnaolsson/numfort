numfort_add_external(fitpack
    fitpack_real64.f
    constraints_tree.f90
)

# This replacement is modern Fortran, unlike the original fixed-form FITPACK.
set_property(SOURCE constraints_tree.f90 APPEND PROPERTY COMPILE_OPTIONS
    "$<$<COMPILE_LANG_AND_ID:Fortran,GNU>:-std=f2018>")
