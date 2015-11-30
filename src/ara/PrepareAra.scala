package ara

import java.io.File
import java.io.PrintWriter

object PrepareAra {

  case class Config(val directory: File = null)

  def main(args: Array[String]) {

    val parser = new scopt.OptionParser[Config]("java -jar ara.jar prepare-ara") {
      opt[File]('d', "directory") required () action { (x, c) => c.copy(directory = x) } text ("Directory with samples.")
    }

    parser.parse(args, Config()) map { config =>

      val dir = config.directory
      val samples = dir.listFiles().filter(_.isDirectory())

      //var count = 0
      
      /*val jpw1 = new PrintWriter(new File(dir.getPath + "/ara-snptyper.sh"))
      jpw1.println("#!/bin/bash")
      jpw1.println("#SBATCH --job-name=ara-snptyper")
      jpw1.println("#SBATCH --workdir=" + dir)
      jpw1.println("#SBATCH --partition=long --qos=long")
      jpw1.println("#SBATCH --output=ara-snptyper.out")
      jpw1.println*/
      
      /*val jpw2 = new PrintWriter(new File(dir.getPath + "/ara-interpreter.sh"))
      jpw2.println("#!/bin/bash")
      jpw2.println("#SBATCH --job-name=ara-interpreter")
      jpw2.println("#SBATCH --workdir=" + dir)
      jpw2.println("#SBATCH --partition=long --qos=long")
      jpw2.println("#SBATCH --output=ara-interpreter.out")
      jpw2.println*/

      samples.foreach { s =>
        val name = s.getName
        val bams = s.listFiles().filter(f => f.getName.endsWith("sorted.bam") && !f.getName.contains("dr-region")).map(_.getName)

        if (!bams.isEmpty) {
          
          /*if (count == 100) {
            jpw1.println("sleep 30m")
          }*/
          
          /** SNP-typer sbatch job */
          /*val pw1 = new PrintWriter(new File(s.getPath + "/" + name + ".ara-snptyper.sh"))
          pw1.println("#!/bin/bash")
          pw1.println("#SBATCH --job-name=ara-snptyper-" + name)
          pw1.println("#SBATCH --workdir=" + s)
          pw1.println("#SBATCH --partition=short --qos=short")
          pw1.println("#SBATCH --mem=3000")
          pw1.println("#SBATCH --output=ara-snptyper.out")
          pw1.println
          pw1.println("java -Xmx3g -jar /tudelft.net/staff-bulk/ewi/insy/DBL/Arlin/ara-development/ara.jar snp-typer -o " + name + " " + bams.mkString(" "))
          pw1.close*/

          /** Ara interpreter sbatch job */
          val pw2 = new PrintWriter(new File(s.getPath + "/" + name + ".ara-interpreter.sh"))
          pw2.println("#!/bin/bash")
          pw2.println("#SBATCH --job-name=ara-interpreter-" + name)
          pw2.println("#SBATCH --workdir=" + s)
          pw2.println("#SBATCH --partition=short --qos=short")
          pw2.println("#SBATCH --mem=1000")
          pw2.println("#SBATCH --output=ara-interpreter.out")
          pw2.println
          pw2.println("java -jar /tudelft.net/staff-bulk/ewi/insy/DBL/Arlin/ara-development/ara.jar interpret -m " + name + ".ara -o " + name)
          
          pw2.close
          
          //jpw1.println("sbatch " + name + "/" + name + ".ara-snptyper.sh")
          //jpw2.println("sbatch " + name + "/" + name + ".ara-interpreter.sh")
          
          //count = count + 1
        }

      }

      //jpw1.close
      //jpw2.close

    }
  }
}