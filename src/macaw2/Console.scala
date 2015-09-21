package macaw2

object Console {
  
  def main(args: Array[String]) ={
    
    if (args.size == 0) listInstructions
    else {
      args(0) match {
        
        case "help" => listInstructions
        case "SNP-typer" => MacawSNPtyper.main(args.drop(1))
        case "interpret-MI" => MacawUtilities.main(args.drop(1))
        case "interpret-DR" => DrugResistances.main(args.drop(1))
        case "mutation-rate-DR" => MutationRate.main(args.drop(1))
        case "vcf2snp-phylip" => Vcf2snpPhylip.main(args.drop(1))
        case "qc-snp-sequences" => QCsnpSequences.main(args.drop(1))
        case "qc-snp-fasta" =>QCsnpFasta.main(args.drop(1))
        case "dr-markers" => DrugResistanceMarkers.main(args.drop(1))
        case _ => listInstructions 
      }
    }
    
    
    
    def listInstructions() {
      println("Usage: java -jar ara.jar [instruction] [instruction options...]")
      println("Instructions:")
      println("\tSNP-typer\t\tDetect presence/absence of SNP markers to detect mixed infection or drug resistance.")
      println("\tinterpret-MI\t\tInterpret MTBC hierarchy clusters from SNP-typer results.")
      println("\tinterpret-DR\t\tInterpret drug resistances from SNP-typer results.")
      println("\tvcf2snp-phylip\t\tWrite phy-file from SNP sequences.")
      println("\tmutation-rate-DR\t\tCount number of mutations within 21 bp window of DR markers given of list of VCFs.")
      println("\tqc-snp-sequences\t\tCheck whether there are snp-sequences in a phy-file with more than 5% of N's in the total sequence.")
      println("\tqc-snp-fasta\t\tCheck length or N-content of multi fasta-file with SNP sequences.")
      println("\tdr-markers\t\tCreate drug resistance and susceptibility markers from a list with TB drug resistances.")
    }    
    
  }
  
}