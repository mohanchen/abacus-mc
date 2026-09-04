#ifndef OCC_MIXER_H
#define OCC_MIXER_H

#include "source_estate/occ_matrix.h"

#include <vector>

class UnitCell;
class Charge_Mixing;

/**
 * @brief Mixing of the DFT+U on-site occupation matrix.
 *
 * Owns the flattened occupation-matrix buffers used by the charge-mixing
 * machinery (PW path) and the plain linear mixing kernel (LCAO path).
 *
 * The flat layout reuses the pot_uterm_pw_index offset table: for nspin=2
 * the buffer is split into [all_up | all_dn] halves; for nspin=1/4 a single
 * block per atom is used. Serialization to/from the nested OccupationMatrix
 * is delegated to OccupationMatrix::{write_to_flat, read_from_flat,
 * write_save_to_flat}.
 *
 * An OccMatMixer instance exists only when mixing is enabled, so its
 * presence doubles as the "mixing on" flag (no mutable workflow switch).
 */
class OccMatMixer
{
  public:
    OccMatMixer() = default;
    ~OccMatMixer() = default;

    /**
     * @brief Allocate the flat buffers and bind the layout table.
     * @param cell          unit cell (borrowed, must outlive this object)
     * @param orbital_corr  per-type correlated-l table (borrowed)
     * @param flat_index    per-atom offset table, i.e. pot_uterm_pw_index (borrowed)
     * @param nspin         spin channels (1, 2 or 4)
     * @param total_size    total flat-buffer size (== pot_uterm_pw.size())
     */
    void init(const UnitCell* cell,
              const std::vector<int>* orbital_corr,
              const std::vector<int>* flat_index,
              int nspin,
              int total_size);

    /**
     * @brief Seed uom_save from an occupation matrix loaded from file.
     *
     * Used when occ_mat_ctrl != 0 (restart from dm_onsite_ini.txt) so that
     * the first mixing step has a meaningful "previous" matrix.
     */
    void seed_save(const OccupationMatrix& occmat);

    /**
     * @brief Begin an SCF iteration: snapshot occ into save buffers.
     *
     * Performs OccupationMatrix::copy_to_save and flattens occ_save into
     * uom_save_.
     */
    void begin_iter(OccupationMatrix& occmat);

    /**
     * @brief Flatten the freshly-computed occupation matrix into uom_.
     */
    void collect(const OccupationMatrix& occmat);

    /**
     * @brief Mix uom_ against uom_save_ via Charge_Mixing (PW path), then
     *        write the mixed result back into the occupation matrix.
     * @param chgmix  charge-mixing object (must have allocate_mixing_uom called)
     */
    void mix(OccupationMatrix& occmat, Charge_Mixing* chgmix);

    /// Total flat-buffer size (== Charge_Mixing::allocate_mixing_uom argument).
    int flat_size() const { return static_cast<int>(uom_.size()); }

  private:
    std::vector<double> uom_;      ///< new / mixed flat occupation matrix
    std::vector<double> uom_save_; ///< previous flat occupation matrix
    const std::vector<int>* index_ = nullptr;        ///< borrowed pot_uterm_pw_index
    const UnitCell* cell_ = nullptr;                 ///< borrowed unit cell
    const std::vector<int>* orbital_corr_ = nullptr; ///< borrowed correlated-l table
    int nspin_ = 0;
};

#endif
