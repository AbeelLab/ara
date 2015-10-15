package ara

import java.io.File
import java.io.PrintWriter
import scala.io.Source

object PrepareDRMapping {

  case class Config(
    val directory: File = null)

  def main(args: Array[String]) {

    val parser = new scopt.OptionParser[Config]("java -jar ara.jar") {
      opt[File]('d', "directory") required () action { (x, c) => c.copy(directory = x) } text ("Directory with samples.")
    }

    parser.parse(args, Config()) map { config =>

      def listFiles(f: Any): List[File] = f match {
        case f: File if (f.isDirectory()) => f.listFiles().toList.flatMap(listFiles(_))
        case f: File if (f.isFile() && f.getName.equals("reduced.vcf")) => List(f)
        case _ => Nil

      }

      val dir = config.directory
      val mainpw = new PrintWriter(new File(dir.getPath + "/dr-region-TU.sh"))
      mainpw.println("#!/bin/bash")
      mainpw.println
      
      val samples = listFiles(dir)
      samples.foreach { rvcf =>
        val sample = rvcf.getParentFile
        val bamfiles = sample.list().filter(_.endsWith(".sorted.bam"))
        val cmd = Source.fromFile(rvcf).getLines.toList(3).dropRight(1).split(" ").drop(8).grouped(2).toList.map{b => 
          val bam = b(1).split("/").last
          b(0) + " " + bam.dropRight(11) + ".dr-region" + bam.takeRight(11)
        }.mkString(" ")
        mainpw.println("sbatch ./" + sample.getName + "/dr-region." + sample.getName + ".sh")
        val pw = new PrintWriter(new File(sample + "/dr-region." + sample.getName + ".sh"))
        pw.println("#!/bin/bash")
        pw.println("#SBATCH --job-name=dr-region_" + sample.getName)
        pw.println("#SBATCH --workdir=" + sample)
        pw.println("#SBATCH --partition=short --qos=short")
        pw.println("#SBATCH --mem=15000")
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
        pw.close
      }

      mainpw.close
      
      println("Finished preparing SLURM scripts for " + dir)
      
    }

  }
}