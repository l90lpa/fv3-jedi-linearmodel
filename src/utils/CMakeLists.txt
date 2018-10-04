list( APPEND utils_files
fv3jedi_lm_const_mod.F90
fv3jedi_lm_kinds_mod.F90
fv3jedi_lm_utils_mod.F90
MAPL_Constants.F90
tapenade/adBuffer.f
tapenade/adStack.c
tapenade/tapenade_iter.F90
)


PREPEND( _p_utils_files       "utils"       ${utils_files} )

set ( utils_src_files
      ${_p_utils_files}
      PARENT_SCOPE
)
