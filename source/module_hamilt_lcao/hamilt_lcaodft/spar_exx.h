#ifdef __EXX
    template<typename Tdata> void cal_HR_exx_sparse(
            const int &current_spin,
            const double &sparse_thr,
            const int (&nmp)[3],
            const std::vector< std::map <int, std::map < std::pair<int, std::array<int,3>>, RI::Tensor<Tdata> > >>& Hexxs);
#endif
