# AutoFE

This code is part of the publication "Data-driven discovery of cardiolipin-selective small molecules by computational active learning", [https://doi.org/10.1039/D2SC00116K](https://doi.org/10.1039/D2SC00116K).
AutoFE is an automated workflow designed to calculate solvation free energies of small molecules in two compared environments. It includes an early-exit strategy if a small-molecule candidate does not show target properties.
The components handling the simulation setup can be run in sequence or independently.

## Installation

Clone the repository and install dependencies:

```bash
git clone https://github.com/Bernadette-Mohr/AutoFE.git
cd AutoFE
pip install -r requirements.txt
```

## Requirements

- Python >= 3.8
- pandas
- scikit-learn
- numpy

(See `requirements.txt` for full dependency list.)

## Contributing

Contributions are welcome! Please see [CONTRIBUTING.md](CONTRIBUTING.md) for details.

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

## Contact

For questions, issues, or feature requests, please open an [issue](https://github.com/Bernadette-Mohr/AutoFE/issues) or contact [Bernadette-Mohr](https://github.com/Bernadette-Mohr).

