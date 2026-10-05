#ifndef DFTU_BASE_IO_H
#define DFTU_BASE_IO_H

#include "source_base/matrix.h"
#include "source_estate/occ_matrix.h"

#include <iosfwd>
#include <string>
#include <vector>

class Plus_U_Base;
class UnitCell;

namespace DFTU_BASE
{

/// nested occupation-matrix type used by DFT+U: occ_mat[iat][l][spin](m0, m1)
using OccMatData = std::vector<std::vector<std::vector<ModuleBase::matrix>>>;

/// Frequency configuration for numbered occupation-matrix output.
///
/// Each output ionic step (out_freq_ion) owns one file dm_onsiteg{#}.txt.
/// Within that file, one section is appended at every electronic step
/// selected by out_freq_elec, scf_nmax or convergence.
struct OccmatOutputCfg
{
    int out_freq_ion;  ///< ionic-step interval; 0 disables numbered files
    int out_freq_elec; ///< electronic-iteration interval recorded inside a numbered file
    int scf_nmax;      ///< maximum number of electronic iterations
};

/// Text format used by write_occup_m().
enum OccmatTextFormat
{
    OCMAT_FMT_LEGACY,   ///< whitespace-token layout parsed by read_occup_m()
    OCMAT_FMT_READABLE  ///< compact human-readable layout for snapshot files
};

/// Storage layout of the real (2l+1)*npol matrix when nspin == 4.
enum OccmatSocLayout
{
    /// 4 contiguous Pauli-component blocks [charge, sigma_x, sigma_y, sigma_z],
    /// each (2l+1)x(2l+1); used by the PW path
    SOC_LAYOUT_PAULI,
    /// real symmetric matrix in spin basis (m + spin-polarization index);
    /// imaginary parts are not stored; used by the LCAO path
    SOC_LAYOUT_SPIN_BASIS_REAL
};

/// Tell whether the given ionic step produces a numbered file.
///
/// @param istep ionic-step index, starting from 0
/// @param cfg frequency configuration
/// @return true when out_freq_ion > 0 and istep is an output ionic step
bool is_ion_step_output_step(int istep, const OccmatOutputCfg& cfg);

/// Tell whether the current electronic step is recorded inside the numbered file.
///
/// @param iter electronic-iteration index, starting from 1
/// @param conv_esolver whether the electronic SCF is converged at this step
/// @param cfg frequency configuration
/// @return true when iter % out_freq_elec == 0, iter == scf_nmax, or converged
bool is_elec_snapshot_trigger(int iter,
                              bool conv_esolver,
                              const OccmatOutputCfg& cfg);

/// Build the name of the per-ionic-step occupation-matrix file.
///
/// @param out_dir output directory (including the trailing separator)
/// @param istep ionic-step index, starting from 0
/// @return dm_onsiteg{istep+1}.txt
std::string gen_ion_step_dm_onsite_filename(const std::string& out_dir, int istep);

/// Append one electronic-step section to the per-ionic-step file.
///
/// The section records the electronic-step index, total energy (converted to
/// eV), total magnetism (Bohr magneton per cell) and convergence status,
/// followed by the occupation matrices and per-atom magnetism. The caller is
/// responsible for truncating the file at the first electronic step of the
/// ionic step (see output()) and for invoking this function only after the
/// total energy and magnetism of the current electronic step are updated.
///
/// @param dftu DFT+U object holding the occupation matrices
/// @param ucell unit cell
/// @param global_out_dir output directory (including the trailing separator)
/// @param nspin number of spin components (1, 2 or 4)
/// @param npol number of polarizations
/// @param istep ionic-step index, starting from 0
/// @param iter electronic-iteration index, starting from 1
/// @param conv_esolver whether the electronic SCF is converged at this step
/// @param etot_ry total energy of the current electronic step, in Ry
/// @param tot_mag total collinear magnetism (Bohr mag/cell)
/// @param tot_mag_nc three non-collinear magnetism components (Bohr mag/cell);
///        may be null when nspin != 4
/// @param cfg frequency configuration
/// @param soc_layout storage layout of the nspin == 4 occupation matrix
void append_ion_step_snapshot(const Plus_U_Base& dftu,
                              const UnitCell& ucell,
                              const std::string& global_out_dir,
                              int nspin,
                              int npol,
                              int istep,
                              int iter,
                              bool conv_esolver,
                              double etot_ry,
                              double tot_mag,
                              const double* tot_mag_nc,
                              const OccmatOutputCfg& cfg,
                              OccmatSocLayout soc_layout);

/// Read the local occupation number matrix from file (rank 0 only).
///
/// The file format matches the output of write_occup_m(). When the file can
/// not be opened, the run quits with an error message that depends on
/// occ_mat_ctrl and init_chg.
void read_occup_m(const UnitCell& ucell,
                  OccupationMatrix& occ,
                  const std::vector<int>& l_channel,
                  const int occ_mat_ctrl,
                  const std::string& fn,
                  const std::string& init_chg,
                  int nspin,
                  int npol);

/// Broadcast the local occupation number matrices from rank 0 to all ranks.
///
/// Implemented in dftu_base_io.cpp (only available in MPI builds).
void local_occup_bcast(const UnitCell& ucell,
                       OccupationMatrix& occ,
                       const std::vector<int>& l_channel,
                       int nspin,
                       int npol);

/// Output DFT+U information (Hubbard U/J, local occupation matrices) to the
/// running log and, when out_chg is set, to disk.
///
/// The file dm_onsite.txt is always overwritten with the latest occupation
/// matrix (used by init_chg=file and NSCF restarts). When cfg.out_freq_ion
/// is positive and istep is an output ionic step, the per-ionic-step file
/// dm_onsiteg{istep+1}.txt is created: at the first electronic step
/// (iter == 1) it is truncated and initialized with a header, then
/// append_ion_step_snapshot() appends one section per recorded electronic
/// step.
///
/// Extracted from Plus_U_Base::output as a free function so that IO logic is
/// decoupled from the Plus_U_Base class. The function only reads the
/// Plus_U_Base state via public accessors; no friend declaration needed.
void output(const Plus_U_Base& dftu,
            const UnitCell& ucell,
            bool out_chg,
            const std::string& global_out_dir,
            int nspin,
            int npol,
            int istep,
            int iter,
            const OccmatOutputCfg& cfg,
            OccmatSocLayout soc_layout);

/// Write local occupation matrices to the given stream.
///
/// Extracted from Plus_U_Base::write_occup_m. When diag is true, eigenvalues
/// and magnetism are also printed; otherwise only raw matrix elements.
/// fmt selects between the legacy token layout (required for dm_onsite.txt
/// because read_occup_m() parses it) and the readable snapshot layout.
/// soc_layout tells how the nspin == 4 storage is arranged (PW Pauli blocks
/// or LCAO real spin-basis matrix).
/// Caller is responsible for opening/closing the stream.
void write_occup_m(const Plus_U_Base& dftu,
                   const UnitCell& ucell,
                   std::ofstream& ofs,
                   bool diag,
                   int nspin,
                   int npol,
                   OccmatTextFormat fmt,
                   OccmatSocLayout soc_layout);

} // namespace DFTU_BASE

#endif
