! Function to find free unit number
INTEGER FUNCTION get_unit_number()
  IMPLICIT NONE

  INTEGER :: lun
  LOGICAL :: file_open

  lun = 9
  file_open = .TRUE.
  DO WHILE ( file_open )
    lun = lun + 1
    INQUIRE( lun, OPENED = file_open )
  END DO
  get_unit_number = lun
  
   RETURN
END FUNCTION get_unit_number 