

void cal_HSR(
		const int &current_spin, 
		const double &sparse_thr, 
		const int (&nmp)[3], 
		hamilt::Hamilt<std::complex<double>>* p_ham);

void cal_HContainer_d(
		const int &current_spin, 
		const double &sparse_threshold, 
		const hamilt::HContainer<double>& hR, 
		std::map<Abfs::Vector3_Order<int>, std::map<size_t, std::map<size_t, double>>>& target);

void cal_HContainer_cd(
		const int &current_spin, 
		const double &sparse_threshold, 
		const hamilt::HContainer<std::complex<double>>& hR, 
		std::map<Abfs::Vector3_Order<int>, 
		std::map<size_t, std::map<size_t, std::complex<double>>>>& target);
