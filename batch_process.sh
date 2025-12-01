set -e

# Usage: ./batch_process.sh /path/to/dir/with/stps
DIR="${1:-.}"

for step in "$DIR"/*.stp; do
    # Skip if no .stp files
    [ -e "$step" ] || continue

    echo "Processing: $step"
    # Uses your argparse main:
    #   brep = BRepMesh(args.step, ...)
    #   ...
    #   viz.start()
    python VTKTool.py --step "$step"
    # Script blocks here until you close the VTK window

    echo "Finished: $step"
done

echo "All STEP files processed."