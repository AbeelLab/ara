package ara

import java.io.File
import java.io.PrintWriter
import scala.io.Source

object ListAraResults {
  case class Config(val directory: File = null, val out: String = null)

  def main(args: Array[String]) {

    val parser = new scopt.OptionParser[Config]("java -jar ara.jar list-results") {
      opt[File]('d', "directory") required () action { (x, c) => c.copy(directory = x) } text ("Data directory")
      opt[String]('o', "output") required () action { (x, c) => c.copy(out = x + ".txt") } text ("Output name.")
    }

    parser.parse(args, Config()) map { config =>
      val study = config.directory.listFiles().filter(f => f.isDirectory())
      val pw = new PrintWriter(new File(config.out))
      pw.println("# Ara results")
      pw.println("# Command: " + args.mkString(" "))
      
      
      study.foreach{study => 
        println("Reading " + study.getName + "...")
        val sampleDirs = study.listFiles().filter(f => f.isDirectory())
        sampleDirs.foreach{ s =>
          val files = s.listFiles().filter(f => f.getName.endsWith(".interpret.ara"))
          if (!files.isEmpty) {
            val line = Source.fromFile(files.head).getLines().toList(5)
            pw.println(study.getName + "\t" + s.getName() + "\t" + line.drop(20))
          }
        }
      }

      pw.close
      
    }
  }
}