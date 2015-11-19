package ara

import java.io.File
import java.io.PrintWriter
import scala.io.Source
import ara.Cluster._

object PrepareBlueJay {

  case class Config(
    val input: File = null,
    //val vcfpaths: File = null,
    val exclude: File = null)

  def main(args: Array[String]) {

    val parser = new scopt.OptionParser[Config]("java -jar ara.jar prepare-bluejay") {
      opt[File]('i', "input") required () action { (x, c) => c.copy(input = x) } text ("Input directory with hierarchical cluster files.")
      opt[File]("exclude") action { (x, c) => c.copy(exclude = x) } text ("File with cluster to exclude.")
    }

    def listFiles(f: Any): List[File] = f match {
      case f: File if (f.isDirectory()) => f.listFiles().toList.flatMap(listFiles(_))
      case f: File if (f.isFile() && f.getName.endsWith("_cluster")) => List(f)
      case _ => Nil
    }

    parser.parse(args, Config()) map { config =>

      val clusterFiles = listFiles(config.input)
      if (clusterFiles.isEmpty) {println("Directory does not contain cluster files. Exiting..."); System.exit(-1)} else {println(clusterFiles.size + " cluster files")}
      
      val exclude = if (config.exclude != null) Source.fromFile(config.exclude).getLines.toList else List.empty[String]

      new File("sbatch-scripts").mkdir()
      new File("lineage-matrices").mkdir()
      new File("associated-snps").mkdir()

      val bjpw = new PrintWriter("bluejay.sh")
      val lpw = new PrintWriter("lineage-matrices.sh")
      
      lpw.println("#!/bin/bash")
      lpw.println("#SBATCH --job-name=lineagematrices")
      lpw.println("#SBATCH --output=slurm-lineagematrices")
      lpw.println("#SBATCH --workdir=/tudelft.net/staff-bulk/ewi/insy/DBL/Arlin/Results/5992taxa/BlueJay/lineage-matrices")
      lpw.println("#SBATCH --partition=short --qos=short")
      lpw.println
      
      clusterFiles.foreach { cFile =>

        val cName = cFile.getName.dropRight(8)

        if (!exclude.contains(cName)) {
          
          val samples = Source.fromFile(cFile).getLines.filterNot(_.startsWith("#")).size
          
          lpw.println("java -jar ../../../../bluejay-development/bluejay.jar lineage-matrix -i ../" + cFile + " -o " + cName + "_lineagematrix -c 1")          
          
          val spw = new PrintWriter("sbatch-scripts/" + cName + "_sbatch")
          spw.println("#!/bin/bash")
          spw.println("#SBATCH --job-name=" + cName)
          spw.println("#SBATCH --workdir=/tudelft.net/staff-bulk/ewi/insy/DBL/Arlin/Results/5992taxa/BlueJay/associated-snps")
          spw.println("#SBATCH --output=slurm_" + cName + ".out")
          spw.println("#SBATCH --partition=short --qos=short")
          spw.println("#SBATCH --mem=11000")
          spw.println
          spw.println
          spw.print("java -jar -Xmx11g ../../../../bluejay-development/bluejay.jar associate --lineage ../lineage-matrices/" + cName + "_lineagematrix --variant ../5992taxa_variantmatrix_transposed.txt -o " + cName + "_ --min-tnr 0.95 --min-tpr 0.95 --min-ppv 0.95 --min-npv 0.95")
          if (cName.hasReference) spw.print(" --absence")
          spw.close

          bjpw.println("sbatch sbatch-scripts/" + cName + "_sbatch")

        }

      }

      bjpw.close
      lpw.close

    }

  }

}