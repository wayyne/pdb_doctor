#!/bin/bash

# Check for required arguments
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <conda|venv>"
    exit 1
fi

# Set the environment tool (conda or venv)
TOOL="$1"
PYTHON_V=3.10.12

# Function to create and set up an environment using conda
setup_conda_env() {
    local env_name="$1"
    local requirements_file="$2"
    
    echo "Setting up Conda environment: $env_name"
    conda create --yes --name "$env_name" python=${PYTHON_V}
    if [ $? -ne 0 ]; then
        echo "Failed to create Conda environment: $env_name"
        exit 1
    fi
    
    echo "Activating environment: $env_name"
    source "$(conda info --base)/etc/profile.d/conda.sh"
    conda activate "$env_name"
    
    echo "Installing requirements from $requirements_file"
    pip install --upgrade pip
    pip install -r "$requirements_file"

    conda install -c conda-forge pdbfixer

    echo "Deactivating Conda environment: $env_name"
    conda deactivate
}

# Function to create and set up an environment using Python venv
setup_venv_env() {
    local env_name="$1"
    local requirements_file="$2"
    
    echo "Setting up Python venv environment: $env_name"
    python3 -m venv "$env_name"
    if [ $? -ne 0 ]; then
        echo "Failed to create Python venv environment: $env_name"
        exit 1
    fi
    
    echo "Activating environment: $env_name"
    source "$env_name/bin/activate"
    
    echo "Installing requirements from $requirements_file"
    pip install --upgrade pip
    pip install -r "$requirements_file"
   
    git clone https://github.com/openmm/pdbfixer.git
    cd pdbfixer && python setup.py install && cd -

    echo "Deactivating Python venv environment: $env_name"
    deactivate
}

# Main script logic
for file in *_env.txt; do
    # Extract environment name from the file name
    env_name="${file%_env.txt}"
    
    # Call the appropriate setup function based on the tool
    case "$TOOL" in
        conda)
            setup_conda_env "$env_name" "$file"
            ;;
        venv)
            setup_venv_env "$env_name" "$file"
            ;;
        *)
            echo "Unsupported tool: $TOOL"
            exit 1
            ;;
    esac
    
    echo "Environment '$env_name' setup complete."
done

echo -e "\n----------------------------------------------\n"
echo "All environments have been set up successfully."
echo "You can activate them and start using them. Cheers!"

