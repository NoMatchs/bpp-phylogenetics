#include <iostream>
#include <memory>
#include <Bpp/Seq/Alphabet/DNA.h>
#include <Bpp/Seq/Io/Fasta.h>
#include <Bpp/Seq/Container/VectorSiteContainer.h>
#include <Bpp/Seq/Container/SiteContainerTools.h> 
#include <Bpp/Phyl/Io/Newick.h>
#include <Bpp/Phyl/Tree/PhyloTree.h>
#include <Bpp/Phyl/Likelihood/ParametrizablePhyloTree.h>
#include <Bpp/Phyl/Model/Nucleotide/GTR.h>
#include <Bpp/Numeric/Prob/GammaDiscreteDistribution.h>
#include <Bpp/Phyl/Likelihood/RateAcrossSitesSubstitutionProcess.h>
#include <Bpp/Phyl/Likelihood/DataFlow/LikelihoodCalculationSingleProcess.h>
#include <Bpp/Phyl/Likelihood/PhyloLikelihoods/SingleProcessPhyloLikelihood.h>

using namespace bpp;

int main() {
    auto dna = std::make_shared<DNA>();
    Fasta fasta;
    auto seqs = fasta.readSequences("real_alignment.fasta", dna);
    if (!seqs) { std::cerr << "读取序列失败" << std::endl; return 1; }
    
    auto rawSites = std::make_shared<VectorSiteContainer>(*seqs);
    std::cout << "原始位点数: " << rawSites->getNumberOfSites() << std::endl;

    // 提取完全位点（返回 unique_ptr）
    auto completeSitesPtr = SiteContainerTools::getCompleteSites(*rawSites);
    
    // 关键修正：用 std::move 转移所有权给 shared_ptr
    auto sites = std::shared_ptr<VectorSiteContainer>(std::move(completeSitesPtr));
    
    std::cout << "完全位点数: " << sites->getNumberOfSites() << std::endl;

    Newick newick;
    auto tree = newick.readPhyloTree("ml_real_output.tree");
    if (!tree) { std::cerr << "读取树失败" << std::endl; return 1; }

    auto model = std::make_shared<GTR>(dna);
    auto rateDist = std::make_shared<GammaDiscreteDistribution>(4, 0.5);
    auto pTree = std::make_shared<ParametrizablePhyloTree>(*tree);

    auto process = std::make_shared<RateAcrossSitesSubstitutionProcess>(model, rateDist, pTree);

    Context context;
    auto likCalc = std::make_shared<LikelihoodCalculationSingleProcess>(context, sites, process);
    auto treeLik = std::make_shared<SingleProcessPhyloLikelihood>(context, likCalc);

    std::cout << "似然值: " << treeLik->getValue() << std::endl;
    std::cout << "done" << std::endl;
    _Exit(0);
}
