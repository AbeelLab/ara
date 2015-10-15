package ara

import java.io.File
import java.io.PrintWriter
import scala.io.Source

object PrepareBlueJay {

  case class Config(
    val input: File = null,
    val output: File = null)

  def main(args: Array[String]) {

    val parser = new scopt.OptionParser[Config]("java -jar ara.jar") {
      opt[File]('i', "input") required () action { (x, c) => c.copy(input = x) } text ("Input directory with hierarchical cluster files.")
      opt[File]('o', "output") required () action { (x, c) => c.copy(output = x) } text ("Output directory to store results.")
    }

    def listFiles(f: Any): List[File] = f match {
      case f: File if (f.isDirectory()) => f.listFiles().toList.flatMap(listFiles(_))
      case f: File if (f.isFile() && f.getName.endsWith("_clusters")) => List(f)
      case _ => Nil
    }

    parser.parse(args, Config()) map { config =>

      val clusterFiles = listFiles(config.input)
      val outputDir = config.output
      
      val pw1 = new PrintWriter(new File("bluejay-linMatrix.sh"))
      pw1.println("#!/bin/bash")
      pw1.println("#SBATCH --job-name=BJ_linMatrix")
      pw1.println("#SBATCH --workdir=" + outputDir)
      pw1.println("#SBATCH --partition=short --qos=short")
      pw1.println("#SBATCH --mem=10000")
      pw1.println
      clusterFiles.foreach { cFile =>
        pw1.println("java -jar ../../../bluejay-3/bluejay.jar lineage-matrix -i " + cFile + " -o " + cFile.getName.dropRight(9) + "_lineagematrix -c 1")
      }      
      pw1.close
      
      val pw2 = new PrintWriter(new File("bluejay-association.sh"))
      pw2.println("#!/bin/bash")
      pw2.println("#SBATCH --job-name=BJ_association")
      pw2.println("#SBATCH --workdir=" + outputDir)
      pw2.println("#SBATCH --partition=long --qos=long")
      pw2.println("#SBATCH --mem=10000")
      pw2.println
      clusterFiles.foreach { cFile =>
        pw2.println("java -jar ../../../bluejay-3/bluejay.jar associate --lineage " + cFile.getName.dropRight(9) + "_lineagematrix --variant 5992taxa_variantmatrix -o " + cFile.getName.dropRight(9) + "_")
      }      
      pw2.close

    }

  }

}