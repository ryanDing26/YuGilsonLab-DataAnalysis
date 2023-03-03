if [ $# -eq 0 ]; then
    echo "Usage: $0 <argument>"
    exit 1
fi

python RMSD.py "/Users/yashravipati/Downloads/PDBBind_processed/$1/$1_docked.pdbqt" "/Users/yashravipati/Downloads/$1.pdb"