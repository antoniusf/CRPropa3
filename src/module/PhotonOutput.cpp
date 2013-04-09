#include "mpc/module/PhotonOutput.h"

#include <iostream>
#include <sstream>

#include "dint/prop_second.h"

// Default cosmological parameters
#define DEFAULT_OMEGA_M 0.3 /**< Default value for _fOmegaM */
#define DEFAULT_OMEGA_LAMBDA 0.7 /**< Default value for _fOmegaLambda */
#define DEFAULT_H_0_KM_S_MPC 71. /**< Default value for _fH0 */

using namespace std;

namespace mpc {

PhotonOutput::PhotonOutput(const string &filename, ref_ptr<MagneticField> field) :
		filename(filename), fout(filename.c_str()), field(field), IRFlag(0), RadioFlag(
				0), H0(DEFAULT_H_0_KM_S_MPC), OmegaLambda(DEFAULT_OMEGA_LAMBDA), OmegaM(
				DEFAULT_OMEGA_M), Zmax(5), Cutcascade_Magfield(0) {
	dataPath = getDataPath("dint");

//TODO: implement ir and radio flags
#if 0
	if (_fpUniverse->Infrared()->Type() == IR_UNIFORM_HIGH)
	IRFlag = 0;
	if (_fpUniverse->Infrared()->Type() == IR_UNIFORM_LOW)
	IRFlag = 1;
	if (_fpUniverse->Infrared()->Type() == IR_UNIFORM_PRIMACK)
	IRFlag = 2;
//  if (lIRFlag != 2) cerr << "Warning : Low/high IR for dint not available otherwise.." << endl;
	if (_fpUniverse->RadioBackground() == "High")
	RadioFlag = 0;
	if (_fpUniverse->RadioBackground() == "Med")
	RadioFlag = 1;
	if (_fpUniverse->RadioBackground() == "Obs")
	RadioFlag = 2;
	if (_fpUniverse->RadioBackground() == "Null")
	RadioFlag = 3;
#endif
}

PhotonOutput::~PhotonOutput() {
}

void PhotonOutput::process(Candidate *candidate) const {
	if (candidate->current.getId() != 22)
		return;

	if (!candidate->hasProperty("Detected"))
		return;

	// Initialize the energy grids for dint
	dCVector energyGrid, energyWidth;
	New_dCVector(&energyGrid, NUM_MAIN_BINS);
	New_dCVector(&energyWidth, NUM_MAIN_BINS);
	SetEnergyBins(MIN_ENERGY_EXP, &energyGrid, &energyWidth);

	// Initialize the spectrum
	Spectrum inputSpectrum;
	NewSpectrum(&inputSpectrum, NUM_MAIN_BINS);
	string spectrum_string;
	if (candidate->getProperty("PhotonSpectrum", spectrum_string)) {
		stringstream spectrum_stream(spectrum_string);
		for (size_t iSp = 0; iSp < NUM_SPECIES; iSp++) {
			for (size_t iB = 0; iB < NUM_MAIN_BINS; iB++) {
				spectrum_stream >> inputSpectrum.spectrum[iSp][iB];
			}
		}
	} else {
		double criticalEnergy = candidate->current.getEnergy()
				/ (eV * ELECTRON_MASS); // units of dint
		int maxBin = (int) ((log10(criticalEnergy * ELECTRON_MASS)
				- MAX_ENERGY_EXP)*BINS_PER_DECADE + NUM_MAIN_BINS);
		inputSpectrum.spectrum[PHOTON][maxBin] = 1.;
	}

	// Initialize the positions
	Vector3d position = candidate->current.getPosition();
	Vector3d initialPosition = candidate->initial.getPosition();
	string initial_string;
	if (candidate->getProperty("PhotonInitialPosition", initial_string)) {
		stringstream initial_stream(initial_string);
		initial_stream >> initialPosition.x >> initialPosition.y
				>> initialPosition.z;
	}

	string initialType = "photon";
	candidate->getProperty("PhotonInitialType", initialType);

	double showerPropDistance = (initialPosition - position).getMag();

	// Initialize the bField
	dCVector bField;
	//TODO: use resolution dependent step size
	int bFieldSteps = showerPropDistance / (0.1 * Mpc);
	New_dCVector(&bField, bFieldSteps);
	Vector3d direction = (initialPosition - position).getUnitVector();
	for (int i = 0; i < bFieldSteps; i++) {
		Vector3d p = initialPosition
				+ direction * showerPropDistance * (i + 0.5) / bFieldSteps;
		Vector3d B = field->getField(p);
		bField.vector[i] = B.perp(direction) / gauss;
	}

	// Initialize output spectrum
	Spectrum outputSpectrum;
	NewSpectrum(&outputSpectrum, NUM_MAIN_BINS);

	prop_second(showerPropDistance / Mpc, &bField, &energyGrid, &energyWidth,
			&inputSpectrum, &outputSpectrum, dataPath, IRFlag, Zmax, RadioFlag,
			H0, OmegaM, OmegaLambda, Cutcascade_Magfield);

#pragma omp critical
	{
		fout << 22 << " " << initialType << " "
				<< candidate->initial.getPosition() / Mpc << " "
				<< initialPosition / Mpc << " " << position / Mpc << " "
				<< candidate->current.getEnergy() << " "
				<< candidate->current.getDirection().x << " "
				<< candidate->current.getDirection().y << " "
				<< candidate->current.getDirection().z << " "
				<< candidate->initial.getEnergy() / Mpc << endl;

		for (int j = 0; j < outputSpectrum.numberOfMainBins; j++) {
			fout << outputSpectrum.spectrum[0][j] << " "; // spectrum: mean number of particles per energy bin
		}
		fout << endl;
	}

	DeleteSpectrum(&outputSpectrum);
	DeleteSpectrum(&inputSpectrum);
	Delete_dCVector(&bField);
	Delete_dCVector(&energyGrid);
	Delete_dCVector(&energyWidth);
}

string PhotonOutput::getDescription() const {
	std::stringstream s;
	s << "PhotonModule: Output file = " << filename;
	return s.str();
}

} // namespace mpc
