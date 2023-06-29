def dock_file(docker_command, protein_filepath, ligand_filepath, center_x, center_y, center_z, size_x, size_y, size_z): # dock the decoy pdbqt to the receptor pdbqt using GWOVina CLI
    os.system(f'{docker_command} --receptor {protein_filepath} --ligand {ligand_filepath}  --center_x  {center_x} --center_y {center_y} --center_z {center_z} --size_x  {size_x} --size_y {size_y}  --size_z {size_z}' \
              ' --exhaustiveness=16 --num_wolves=16 --num_modes=1 --energy_range=3')


def dock_all_decoys(decoy_pdbqts, pdbqt_files, docker_command, docked_decoys, padding): # batch dock all decoy.pdbqt files in 'decoy_pdbqts' folder

    # iterate over pdbqt decoy files and dock
    decoy_pdbqt_files = [f'{decoy_pdbqts}{file}' for file in os.listdir(decoy_pdbqts)]
    with tqdm(total=len(decoy_pdbqt_files)) as pbar:
        for filepath, filename in zip(decoy_pdbqt_files, os.listdir(decoy_pdbqts)):

            # define file variables
            pdb_code = filename.split('_')[0]
            foldername = filename.split('.')[0]
            receptor_file = f'{pdbqt_files}{pdb_code}/{pdb_code}_receptor.pdbqt'
            example_crystal_ligand = f'{pdbqt_files}{pdb_code}/{pdb_code}_ligand.pdbqt'

            # dock the decoy to the active crystal receptor.pdbqt
            try:
                dock_file(docker_command, receptor_file, filepath, *get_coordinates(example_crystal_ligand, padding))

                # make destination folder and transfer AutoDockTools output decoy.pdbqt to destination folder
                os.mkdir(f'{docked_decoys}{foldername}')
                shutil.copyfile(f'{decoy_pdbqts}{foldername}_out.pdbqt',f'{docked_decoys}{foldername}/{foldername}_ligand.pdbqt')
                shutil.copyfile(receptor_file, f'{docked_decoys}{foldername}/{foldername}_receptor.pdbqt')
                os.remove(f'{decoy_pdbqts}{foldername}_out.pdbqt')
            except:
                print(f'ERROR: Could not dock {filename}')
            pbar.update(1)


def parse_args(args): # parse CLI user inputs

    docker_command = args[args.index('-dock') + 1]

    ref_df = pd.read_csv(args[args.index('-ref') + 1])

    pdbqt_files = args[args.index('-pdbqts') + 1]

    decoys = args[args.index('-decoys') + 1]

    docked_decoys = args[args.index('-des') + 1]

    decoy_files = [f'{decoys}{file}' for file in os.listdir(decoys)]

    mgl_tools_path = args[args.index('-mgl') + 1]

    padding = int(args[args.index('-pad') + 1])

    return docker_command, ref_df, pdbqt_files, decoys, docked_decoys, decoy_files, mgl_tools_path, padding


def main(): # run script using CLI

    docker_command, ref_df, pdbqt_files, decoys, docked_decoys, decoy_files, mgl_tools_path, padding = parse_args(sys.argv)

    # make output folders if none exist
    decoy_pdbs = os.path.join(os.path.dirname(__file__), 'decoy_pdbs','')
    decoy_pdbqts = os.path.join(os.path.dirname(__file__), 'decoy_pdbqts','')
    if not os.path.isdir(decoy_pdbs):
        os.mkdir(decoy_pdbs)
    if not os.path.isdir(decoy_pdbqts):
        os.mkdir(decoy_pdbqts)

    # first convert decoys from smiles to pdbqt
    make_pdbs_from_smiles(ref_df, decoy_files, decoy_pdbs)
    make_pdbqts_from_pdbs(decoy_pdbs, decoy_pdbqts, mgl_tools_path)

    # dock the decoy pdbqts to the respetive crystal receptor.pdbqt
    dock_all_decoys(decoy_pdbqts, pdbqt_files, docker_command, docked_decoys, padding)

if __name__ == '__main__':
    main()