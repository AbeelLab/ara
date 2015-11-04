package ara

import java.io.File
import java.io.PrintWriter
import scala.io.Source

object PrepareBlueJay {

  case class Config(
    val input: File = null,
    //val vcfpaths: File = null,
    val exclude: File = null)

  def main(args: Array[String]) {

    val parser = new scopt.OptionParser[Config]("java -jar ara.jar prepare-bluejay") {
      opt[File]('i', "input") required () action { (x, c) => c.copy(input = x) } text ("Input directory with hierarchical cluster files.")
      //opt[File]("vcf-paths") required () action { (x, c) => c.copy(vcfpaths = x) } text ("File with all VCF-paths")
      opt[File]("exclude") action { (x, c) => c.copy(exclude = x) } text ("File with clusters to exclude.")
    }

    def listFiles(f: Any): List[File] = f match {
      case f: File if (f.isDirectory()) => f.listFiles().toList.flatMap(listFiles(_))
      case f: File if (f.isFile() && f.getName.endsWith("_clusters")) => List(f)
      case _ => Nil
    }

    parser.parse(args, Config()) map { config =>

      val clusterFiles = listFiles(config.input)
      //val vcfpaths = Source.fromFile(config.vcfpaths).getLines.map(line => (line.split("/")(8) -> line)).toMap
      val exclude = if (config.exclude != null) Source.fromFile(config.exclude).getLines.toList else List.empty[String]

      //new File("vcfpathsfiles").mkdir()
      //new File("sbatch-scripts").mkdir()
      //new File("variant-matrices").mkdir()

      val bjpw = new PrintWriter("bluejay.sh")

      clusterFiles.foreach { cFile =>

        val cName = cFile.getName.dropRight(8)

        if (!exclude.contains(cName)) {
          val samples = Source.fromFile(cFile).getLines.filterNot(_.startsWith("#")).size
          
          /**
           * val vcfs = Source.fromFile(cFile).getLines.filterNot(_.startsWith("#")).filterNot(_.startsWith("reference")).map(line => vcfpaths(line.split("\t")(0))).toList
           * val vcfpw = new PrintWriter("vcfpathsfiles/" + cName + "vcfpaths")
           * vcfs.foreach(vcfpw.println)
           * vcfpw.close
           */

          val spw = new PrintWriter("sbatch-scripts/" + cName + "sbatch")
          spw.println("#!/bin/bash")
          spw.println("#SBATCH --job-name=" + cName)
          spw.println("#SBATCH --workdir=/tudelft.net/staff-bulk/ewi/insy/DBL/Arlin/Results/6100taxa/BlueJay_5992taxa")
          spw.println("#SBATCH --output=slurm_" + cName + ".out")
          if (samples > 1500) {
            spw.println("#SBATCH --partition=long --qos=long")
          } else {
            spw.println("#SBATCH --partition=short --qos=short")
          }

          val mem = (samples / 200).toInt + 5
          spw.println("#SBATCH --mem=" + mem * 1000)
          spw.println
          spw.println("java -jar ../../../bluejay-development/bluejay.jar variant-matrix -v vcfpathsfiles/" + cFile.getName.dropRight(8) + "vcfpaths -o variant-matrices/" + cName + " --snps-only")
          spw.println("wait")
          spw.println("java -jar -Xmx" + mem + "g ../../../bluejay-development/bluejay.jar associate --lineage lineagematrices/" + cFile.getName.dropRight(8) + "lineagematrix --variant variant-matrices/" + cName + "variantmatrix.txt -o " + cName + "--min-tnr 0.95 --min-tpr 0.95 --min-ppv 0.95 --min-npv 0.95")
          spw.close

          bjpw.println("sbatch sbatch-scripts/" + cName + "sbatch")

        }

      }

      bjpw.close

    }

  }

}