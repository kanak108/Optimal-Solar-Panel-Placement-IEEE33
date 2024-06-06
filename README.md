# Optimal Solar Panel Placement in IEEE 33 Bus System

This project aims to determine the optimal placement of a 1 MW solar panel in the IEEE 33 bus system to minimize power losses. The project includes MATLAB scripts to perform load flow analysis and calculate power losses before and after the addition of the solar panel.

## Table of Contents
- [Introduction](#introduction)
- [Requirements](#requirements)
- [Usage](#usage)
- [Results](#results)
- [Contributing](#contributing)
- [License](#license)

## Introduction
The goal of this project is to find the best bus for placing a 2.5 MW solar panel in the IEEE 33 bus system. By optimizing the placement, we aim to reduce power losses and improve the overall efficiency of the system.

## Requirements
- MATLAB R2018a or later
- Basic knowledge of power systems and load flow analysis

## Usage
1. Clone this repository:
    ```bash
    git clone https://github.com/your-username/Optimal-Solar-Panel-Placement-IEEE33.git
    cd Optimal-Solar-Panel-Placement-IEEE33
    ```

2. Load the bus data and line data by running the provided `Z_bus.m` and `load_data.m` files.

3. Run the main MATLAB script `optimal_placement.m` to find the best bus for solar panel placement:
    ```matlab
    optimal_placement
    ```

4. The script will output:
    - The best bus for solar panel placement.
    - Total real power loss before and after the solar panel addition.
    - Percentage reduction in losses.
    - Voltage profiles before and after the solar panel addition.

## Results
- **Best Bus for Solar Panel Placement:** The bus number with the minimum power loss.
- **Total Real Power Loss Before Addition:** The power loss in the system without the solar panel.
- **Total Real Power Loss After Addition:** The power loss in the system with the solar panel at the best bus.
- **Percentage Reduction in Losses:** The percentage reduction in power loss due to the optimal placement of the solar panel.

## Contributing
Contributions are welcome! If you have any ideas or improvements, feel free to open an issue or submit a pull request.

1. Fork the repository.
2. Create a new branch:
    ```bash
    git checkout -b feature-branch
    ```
3. Commit your changes:
    ```bash
    git commit -m 'Add some feature'
    ```
4. Push to the branch:
    ```bash
    git push origin feature-branch
    ```
5. Open a pull request.

## License
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

