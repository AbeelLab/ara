package ara

import scala.io.Source
import java.io.File
import java.io.PrintWriter

object PreparePilon {
  
  case class Config(
    val input: File = null
  )
  
  def main(args: Array[String]) {
    
    val parser = new scopt.OptionParser[Config]("java -jar ara.jar") {
      opt[File]('s', "samples") required() action { (x, c) => c.copy(input = x) } text ("Input file with samples names.")
    }
    
    parser.parse(args, Config()) map { config =>
      val samples = Source.fromFile(config.input).getLines.filterNot(_.startsWith("#"))
      
      
      samples.foreach{ s => 
        println("dos2unix ./" + s + "/pilon." + s + ".sh")
        println("sbatch ./" + s + "/pilon." + s + ".sh")
        val pw = new PrintWriter("pilon." + s + ".sh" )
        pw.println("#!/bin/bash")
        pw.println("#SBATCH --mem=5000")
        pw.println("#SBATCH --partition=short --qos=short")
        pw.println("#SBATCH -J pilon." + s)
        pw.println("#SBATCH -o pilon.%j.log")
        pw.println("#SBATCH -e pilon.%j.err")
        pw.println("#SBATCH --workdir=/tudelft.net/staff-bulk/ewi/insy/tbarc/data/Lit.Comas2013/" + s)
        pw.println("java -Xmx4g -jar ../../../bin/pilon-1.11.jar --output pilon --vcf --changes --fix all,breaks --genome /tudelft.net/staff-bulk/ewi/insy/tbarc/team/arlin/MT_H37RV_BRD_V5.ref.fasta --unpaired " + s + ".sorted.bam")
        pw.println("wait")
        pw.println("java -Xmx1g -jar ../../../team/arlin/reducedvcf-12/reducevcf.jar -i pilon.vcf -o reduced.vcf -k")
        pw.println("wait")
        pw.println("rm pilon.vcf")
        pw.close
      }
      
      
    }
    
  }
}