def ligand_pose_generator(params, lower_index, upper_index):

    """
    Function: Generate a list of receptor and ligand arguments
              between specific indices to pass to a function to be scored.
              This is necessary because if we do this at the start, if
              many ligands are being scored we cannot store them all in 
              memory and the script crashes

    Inputs:   Command line arguments, the lower and upper indexes
              of ligand poses to score

    Outputs:  List of tuples to score in the form
              [(receptor filename, ligand filename, (ligand pose number, ligand pdbqt block))]
    """

    # set up list to populate
    requested_receptor_ligand_args = list()

    # track the pose index
    pose_index = 0

    # then for each ligand
    for ligand_index, ligand_filepath in enumerate(params.ligand):

         
         # don't waste any time if we're already over the upper index
        if pose_index > upper_index:
            break

        # load the poses from the current ligand being considered
        pdbqt_pose_blocks = list()
        lig_text = open(ligand_filepath, 'r').read()
        lig_poses = lig_text.split('MODEL')
        for pose in lig_poses:
            lines = pose.split('\n')
            clean_lines = [line for line in lines if not line.strip().lstrip().isnumeric() and 'ENDMDL' not in line]
            if len(clean_lines) < 3:
                pass
            else:
                pose = '\n'.join(clean_lines)
                pdbqt_pose_blocks.append(pose)

                # stop if we only want one pose
                if params.pose_1:
                    break

        # make a tuple with pdbqt block and pose name
        poses = [(f'_pose_{pdbqt_pose_blocks.index(pose) + 1}', pose) for pose in pdbqt_pose_blocks]

        # for each pose
        for pose in poses:
            
            # if the pose is one we want
            if lower_index <= pose_index < upper_index:
            
                # add it to the receptor ligand arguments
                receptor_ligand_args = (params.receptor[ligand_index], ligand_filepath, pose)

                requested_receptor_ligand_args.append(receptor_ligand_args)
    
            # update the pose index
            pose_index += 1

    # return the requested poses between the specified indexes
    return requested_receptor_ligand_args
