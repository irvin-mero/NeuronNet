#include "network.h"
#include "random.h"

void Network::resize(const size_t &n, double inhib) {
    size_t old = size();
    neurons.resize(n);
    if (n <= old) return;
    size_t nfs(inhib*(n-old)+.5);
    set_default_params({{"FS", nfs}}, old);
}

void Network::set_default_params(const std::map<std::string, size_t> &types,
                                 const size_t start) {
    size_t k(0), ssize(size()-start), kmax(0);
    std::vector<double> noise(ssize);
    _RNG->uniform_double(noise);
    for (auto I : types) 
        if (Neuron::type_exists(I.first)) 
            for (kmax+=I.second; k<kmax && k<ssize; k++) 
                neurons[start+k].set_default_params(I.first, noise[k]);
    for (; k<ssize; k++) neurons[start+k].set_default_params("RS", noise[k]);
}

void Network::set_types_params(const std::vector<std::string> &_types,
                               const std::vector<NeuronParams> &_par,
                               const size_t start) {
    for (size_t k=0; k<_par.size(); k++) {
        neurons[start+k].set_type(_types[k]);
        neurons[start+k].set_params(_par[k]);
    }
}

void Network::set_values(const std::vector<double> &_poten, const size_t start) {
    for (size_t k=0; k<_poten.size(); k++) 
        neurons[start+k].potential(_poten[k]);
}

bool Network::add_link(const size_t &a, const size_t &b, double str) {
    if (a==b || a>=size() || b>=size() || str<1e-6) return false;
    if (links.count({a,b})) return false;
    if (neurons[b].is_inhibitory()) str *= -2.0;
    links.insert({{a,b}, str});
    return true;
}

size_t Network::random_connect(const double &mean_deg, const double &mean_streng) {
    links.clear();
    std::vector<int> degrees(size());
    _RNG->poisson(degrees, mean_deg);
    size_t num_links = 0;
    std::vector<size_t> nodeidx(size());
    std::iota(nodeidx.begin(), nodeidx.end(), 0);
    for (size_t node=0; node<size(); node++) {
        _RNG->shuffle(nodeidx);
        std::vector<double> strength(degrees[node]);
        _RNG->uniform_double(strength, 1e-6, 2*mean_streng);
        int nl = 0;
        for (size_t nn=0; nn<size() && nl<degrees[node]; nn++)
            if (add_link(node, nodeidx[nn], strength[nl])) nl++;
        num_links += nl;
    }
    return num_links;
}

std::pair<size_t, double> Network::degree(const size_t& rec_neuro) const{
	std::vector<std::pair<size_t, double> > connected_neuros(neighbors(rec_neuro));
	std::pair<size_t,double> res;
	res.first = connected_neuros.size();
	double sum_int(0.0);
	for(const auto& neuro: connected_neuros){
		sum_int+=neuro.second;
	}
	res.second = sum_int;
	return res;
}

std::vector<std::pair<size_t, double> > Network::neighbors(const size_t& ind_recev_neuro) const{
	std:: vector<std::pair<size_t,double>> resul;
	for(linkmap::const_iterator i=links.lower_bound({ind_recev_neuro,0}); i!=links.cend() and (i->first).first == ind_recev_neuro;++i){
		//vu que le deuxieme sur link est le sending et premier receiver on voit que le premier
			resul.push_back(std::pair<size_t,double>(i->first.second,i->second));
	}
	return resul;
	
	
}

std::vector<double> Network::potentials() const {
    std::vector<double> vals;
    for (size_t nn=0; nn<size(); nn++)
        vals.push_back(neurons[nn].potential());
    return vals;
}

std::vector<double> Network::recoveries() const {
    std::vector<double> vals;
    for (size_t nn=0; nn<size(); nn++)
        vals.push_back(neurons[nn].recovery());
    return vals;
}
std::set<size_t> Network::step(const std::vector<double>& thalamic_input){
	std::set<size_t> fir;
	//on prend en compte que ces qui sont firing on garde leur indice pour la suite
	for(size_t i(0);i<neurons.size(); ++i){
		if(neurons[i].firing()){
			fir.insert(i);
			neurons[i].reset();
		}
	}
	
	for(size_t i(0);i<neurons.size(); ++i){
		double to_rem(0.0);
		double to_add(0.0);
		std::vector<std::pair<size_t, double> > connected(neighbors(i));
		//on regarde en premier s'il y a des neurones firing et la neurone i a des neurones connectées,sinon on passe au suivant
		if (not (fir.empty() or connected.empty())){
			for(const auto& connex:connected){
				if(fir.count(connex.first)==1){
					if(neurons[connex.first].is_inhibitory()){
						to_rem+=connex.second;
					}else{ 
						to_add+=connex.second;
					}
				}
			}
		}
		if(neurons[i].is_inhibitory()){
			neurons[i].input(0.4*thalamic_input[i]+0.5*to_add+to_rem);
		} else { neurons[i].input(thalamic_input[i]+0.5*to_add+to_rem);
		}
	
		neurons[i].step();
	}
	
	return fir;
}
void Network::print_params(std::ostream *_out) {
    (*_out) << "Type\ta\tb\tc\td\tInhibitory\tdegree\tvalence" << std::endl;
    for (size_t nn=0; nn<size(); nn++) {
        std::pair<size_t, double> dI = degree(nn);
        (*_out) << neurons[nn].formatted_params() 
                << '\t' << dI.first << '\t' << dI.second
                << std::endl;
    }
}

void Network::print_head(const std::map<std::string, size_t> &_nt, 
                         std::ostream *_out) {
    size_t total = 0;
    for (auto It : _nt) {
        total += It.second;
        for (auto In : neurons)
            if (In.is_type(It.first)) {
                (*_out) << '\t' << It.first << ".v"
                        << '\t' << It.first << ".u"
                        << '\t' << It.first << ".I";
                break;
            }
    }
    if (total<size())
        for (auto In : neurons) 
            if (In.is_type("RS")) {
                (*_out) << '\t' << "RS.v" << '\t' << "RS.u" << '\t' << "RS.I";
                break;
            }
    (*_out) << std::endl;
}

void Network::print_traj(const int time, const std::map<std::string, size_t> &_nt, 
                         std::ostream *_out) {
    (*_out)  << time;
    size_t total = 0;
    for (auto It : _nt) {
        total += It.second;
        for (auto In : neurons) 
            if (In.is_type(It.first)) {
                (*_out) << '\t' << In.formatted_values();
                break;
            }
    }
    if (total<size())
        for (auto In : neurons) 
            if (In.is_type("RS")) {
                (*_out) << '\t' << In.formatted_values();
                break;
            }
    (*_out) << std::endl;
}
