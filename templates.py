def get_header():
    return """'tmp
'File created: Sunday September 01 2019 - 11:36 AM
'SHOPBOT FILE IN INCHES
IF %(25)=1 THEN GOTO UNIT_ERROR	'check to see software is set to standard
C#,90				 	'Lookup offset values
'
'Turning router ON
SO,1,1
PAUSE 2
'
'
'Toolpath Name = Profile 2
'Tool Name   = End Mill (0.25 inches) Carbide
"""

def get_footer(end_x = 0.0, end_y = 0.0):
    return """'
'Turning router OFF
SO,1,0
J2,{0:.6f},{1:.6f}
END
UNIT_ERROR:
C#,91					'Run file explaining unit error
END
""".format(end_x, end_y)