################################################################################
# musica-rad-xfer library

set(RAD_XFER_SRC
    rad_abs_xsect.F90
    rad_trans.F90
    tuv_radiation_transfer.F90
    tuv_subs.F90
   )

add_library(musica-rad-xfer-lib ${RAD_XFER_SRC})

target_link_libraries(musica-rad-xfer-lib ${NETCDF_LIBS})

set_target_properties(musica-rad-xfer-lib PROPERTIES OUTPUT_NAME musica-rad-xfer)

################################################################################
