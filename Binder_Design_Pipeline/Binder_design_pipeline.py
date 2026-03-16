# Install dependencies (skip if already done)
import os

# Set environment variables
os.environ['CCD_MIRROR_PATH'] = ''
os.environ['PDB_MIRROR_PATH'] = ''

if not os.path.isfile("FOUNDRY_READY"):
    print("Installing rc-foundry...")

    # Uninstall torchvision first to avoid operator conflicts
    os.system("pip uninstall -y torchvision")

    # Install rc-foundry
    os.system("pip install -q 'rc-foundry[all]'")

    # Mark as ready
    os.system("touch FOUNDRY_READY")

    print("Done!")
else:
    print("rc-foundry already installed.")