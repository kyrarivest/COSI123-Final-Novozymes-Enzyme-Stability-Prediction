# % $Id: residue_rmsd.tcl,v 1.4 2006/03/06 23:56:46 timisgro Exp $

proc rmsd_residue_over_time {{mol top} res} {
    
    # use frame 0 for the reference
    set reference [atomselect $mol "protein" frame "last"]
    # the frame being compared
    set compare [atomselect $mol "protein"]
    #make a selection with all atoms
    set all [atomselect top all]
    #get the number of frames
    set num_steps [molinfo $mol get numframes]
    #open file for writing
    set fil [open residue_rmsd_sasa_last.dat w]
    
    foreach r $res {
	set rmsd($r) 0
    }
    
    #loop over all frames in the trajectory
    for {set frame 0} {$frame < $num_steps} {incr frame} {
	puts "Calculating rmsd for frame $frame ..."
	# get the correct frame
	$compare frame $frame
        $all frame $frame
	# compute the transformation
	set trans_mat [measure fit $compare $reference]
	# do the alignment
	$all move $trans_mat
	
	# compute the RMSD
	#loop through all residues
	foreach r $res {
	    set ref [atomselect $mol "protein and resid $r and noh" frame "last"]
	    set comp [atomselect $mol "protein and resid $r and noh" frame $frame]
	    set rmsd($r) [expr $rmsd($r) + [measure rmsd $comp $ref]]
	    $comp delete
	    $ref delete
	}
    }
    
    set ave 0
	foreach r $res {
	    set rmsd($r) [expr $rmsd($r)/$num_steps]
	    # SASA using first frame:
	    set selprot [atomselect top "protein" frame 0]
	    set selres [atomselect top "protein and resid $r" frame 0]
	    set sasafirst [measure sasa 1.4 $selprot -restrict $selres]
	    # SASA using last frame:
	    set selprot [atomselect top "protein" frame "last"]
	    set selres [atomselect top "protein and resid $r" frame "last"]
	    set sasalast [measure sasa 1.4 $selprot -restrict $selres]
	    # print the RMSD
	    puts "RMSD of residue $r is $rmsd($r)"
	    puts $fil " $r \t $rmsd($r) \t $sasafirst \t $sasalast"
	    set res_b [atomselect $mol "resid $r"] 
            $res_b set user $rmsd($r)
            $res_b delete
	    set ave [expr $ave + $rmsd($r)]
	}
	
    set ave [expr $ave/[llength $res]]
    puts " Average rmsd per residue:   $ave"
    close $fil
}

