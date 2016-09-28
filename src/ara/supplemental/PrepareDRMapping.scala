package ara.supplemental

import java.io.File
import java.io.PrintWriter
import scala.io.Source

object PrepareDRMapping {

  case class Config(
    val directory: File = null, samplesToSkip: File = null)

  def main(args: Array[String]) {

    val parser = new scopt.OptionParser[Config]("java -jar ara.jar prepare-drmapping") {
      opt[File]('d', "directory") required () action { (x, c) => c.copy(directory = x) } text ("Directory with samples.")
      opt[File]("skip-samples") action { (x, c) => c.copy(samplesToSkip = x) } text ("File with samples to skip.")
    }

    parser.parse(args, Config()) map { config =>

      /*def listReducedVCF(f: Any): List[File] = f match {
        case f: File if (f.isDirectory()) => f.listFiles().toList.flatMap(listReducedVCF(_))
        case f: File if (f.isFile() && f.getName.equals("reduced.vcf")) => List(f)
        case _ => Nil
      }*/

      def listDRVCF(f: Any): List[File] = f match {
        case f: File if (f.isDirectory()) => f.listFiles().toList.flatMap(listDRVCF(_))
        case f: File if (f.isFile() && f.getName.equals("pilon.dr-region.vcf")) => List(f)
        case _ => Nil
      }
      
      val dir = config.directory
      val mainpw = new PrintWriter(new File(dir.getPath + "/dr-region-TU2.sh"))
      mainpw.println("#!/bin/bash")
        mainpw.println("#SBATCH --job-name=readDRvcf")
        mainpw.println("#SBATCH --workdir=" + dir.getPath)
        mainpw.println("#SBATCH --partition=long --qos=long")
        mainpw.println("#SBATCH --mem=10000")
        mainpw.println("#SBATCH --output=dr-region-TU2.out")
        mainpw.println
      mainpw.println
            
      val samples = if (config.samplesToSkip != null) listDRVCF(dir).filterNot(s => Source.fromFile(config.samplesToSkip).getLines.contains(s.getParentFile.getName)) else listDRVCF(dir)
      var count = 0
      samples.foreach { drvcf =>
        count = count + 1
        val sample = drvcf.getParentFile        
        mainpw.println("java -jar /tudelft.net/staff-bulk/ewi/insy/DBL/Arlin/ara-development/ara.jar dr-snps --vcf " + drvcf + " --fasta ../team/arlin/MT_H37RV_BRD_V5.ref.fasta -- gff ../team/arlin/MT_H37RV_BRD_V5.ref.gff -o " + sample.getName)
        if (count % 100 == 0) mainpw.println("sleep 10s")
        /*val bamfiles = sample.list().filter(_.endsWith(".sorted.bam"))
        val cmd = Source.fromFile(drvcf).getLines.toList(3).dropRight(1).split(" ").drop(8).grouped(2).toList.map{b => 
          val bam = b(1).split("/").last
          b(0) + " " + bam.dropRight(11) + ".dr-region" + bam.takeRight(11)
        }.mkString(" ")
        mainpw.println("sbatch ./" + sample.getName + "/dr-region." + sample.getName + ".sh")
        //if (count % 200 == 0) mainpw.println("sleep 90m")
        val pw = new PrintWriter(new File(sample + "/dr-region." + sample.getName + ".sh"))
        pw.println("#!/bin/bash")
        pw.println("#SBATCH --job-name=dr-region_" + sample.getName)
        pw.println("#SBATCH --workdir=" + sample)
        pw.println("#SBATCH --partition=long --qos=long")
        pw.println("#SBATCH --mem=" + (if (bamfiles.size > 2) bamfiles.size * 2000 else 4000))
        pw.println("#SBATCH --output=pilon.dr-region.out")
        pw.println
        pw.println("# Map bam-files")
        bamfiles.foreach { bam =>
          pw.println("samtools bamshuf -uO " + bam + " tmp" + bam.dropRight(11) + " | samtools bam2fq - | gzip - | bwa mem -p -v 0 ../../../team/arlin/MT_H37RV_BRD_V5.ref.dr-region.fasta/MT_H37RV_BRD_V5.ref.dr-region.fasta - | samtools view -b - | samtools sort - " + bam.dropRight(11) + ".dr-region.sorted")
        }
        pw.println("wait")
        pw.println
        pw.println("# Index bam-files")
        bamfiles.foreach { bam =>
          pw.println("samtools index " + bam.dropRight(11) + ".dr-region.sorted.bam")
        }
        pw.println("wait")
        pw.println
        pw.println("# Call variants")
        pw.println("java -Xmx4g -jar ../../../bin/pilon-1.11.jar --output pilon.dr-region --vcf --changes --fix all,breaks --genome /tudelft.net/staff-bulk/ewi/insy/tbarc/team/arlin/MT_H37RV_BRD_V5.ref.dr-region.fasta/MT_H37RV_BRD_V5.ref.dr-region.fasta " + cmd)
        pw.close*/
                
      }

      mainpw.close
      
      println("Finished preparing SLURM script for " + dir)
      
    }

  }
}