MODULE String_Utility

  IMPLICIT NONE


  PRIVATE
  PUBLIC :: StrUpCase
  PUBLIC :: StrLowCase
  PUBLIC :: StrCompress

  CHARACTER( * ), PRIVATE, PARAMETER :: LOWER_CASE = 'abcdefghijklmnopqrstuvwxyz'
  CHARACTER( * ), PRIVATE, PARAMETER :: UPPER_CASE = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ' 


CONTAINS


  FUNCTION StrUpCase ( Input_String ) RESULT ( Output_String )

    CHARACTER( * ), INTENT( IN )     :: Input_String
    CHARACTER( LEN( Input_String ) ) :: Output_String

    INTEGER :: i, n


    Output_String = Input_String

    DO i = 1, LEN( Output_String )

      n = INDEX( LOWER_CASE, Output_String( i:i ) )

      IF ( n /= 0 ) Output_String( i:i ) = UPPER_CASE( n:n )

    END DO

  END FUNCTION StrUpCase


  FUNCTION StrLowCase ( Input_String ) RESULT ( Output_String )

    CHARACTER( * ), INTENT( IN )     :: Input_String
    CHARACTER( LEN( Input_String ) ) :: Output_String

    INTEGER :: i, n


    Output_String = Input_String

    DO i = 1, LEN( Output_String )

      n = INDEX( UPPER_CASE, Output_String( i:i ) )

      IF ( n /= 0 ) Output_String( i:i ) = LOWER_CASE( n:n )

    END DO

  END FUNCTION StrLowCase


  FUNCTION StrCompress( Input_String, n ) RESULT ( Output_String )

    CHARACTER( * ),    INTENT( IN )  :: Input_String
    INTEGER, OPTIONAL, INTENT( OUT ) :: n

    CHARACTER( LEN( Input_String ) ) :: Output_String

    INTEGER, PARAMETER :: IACHAR_SPACE = 32
    INTEGER, PARAMETER :: IACHAR_TAB   = 9

    INTEGER :: i, j
    INTEGER :: IACHAR_Character

    Output_String = ' '

    j = 0

    DO i = 1, LEN( Input_String )

      IACHAR_Character = IACHAR( Input_String( i:i ) )

      IF ( IACHAR_Character /= IACHAR_SPACE .AND. &
           IACHAR_Character /= IACHAR_TAB         ) THEN
        j = j + 1
        Output_String( j:j ) = Input_String( i:i )
      END IF

    END DO

    IF ( PRESENT( n ) ) n = j

  END FUNCTION StrCompress

END MODULE String_Utility
