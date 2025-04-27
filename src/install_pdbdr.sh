echo "This script will install a Python package into your site-packages."
echo "WARNING: If you are not using a virtual environment (e.g., venv, conda),"
echo "this installation will occur at the global level."
echo ""
read -p "Do you want to proceed? (y/N): " user_input

# Check user input
if [[ "$user_input" =~ ^[Yy]$ ]]; then

    echo "Proceeding with installation..."
    python3 -m pip install -e .

    if [[ $? -eq 0 ]]; then
        echo "Package installed successfully."
    else
        echo "An error occurred during installation."
    fi

else
    echo "Installation canceled by the user."
fi

