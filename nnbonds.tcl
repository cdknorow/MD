set f "conn.dat"
set fset "conn_who.dat"
set fp [open $f "w"]
set fpset [open $fset "w"]
set numframes [molinfo top get numframes]


for {set i 0} {$i < $numframes} {incr i} {
	set G [atomselect top "name G" frame $i]
	set C [atomselect top "name C" frame $i]
	set result [measure contacts 1.0 $C $G]
	set asize [llength [lindex $result 0]]
	puts $fp $asize
	puts $fpset $result
    }
close $fp
close $fpset

