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
      pw.println("#Study\tsampleID\tPredicted lineage(s) (frequency, read depth)")
      
      study.foreach{study => 
        println("Reading " + study.getName + "...")
        val sampleDirs = study.listFiles().filter(f => f.isDirectory())
        sampleDirs.foreach{ s =>
          val files = s.listFiles().filter(f => f.getName.endsWith(".interpret.ara"))
          if (!files.isEmpty) {
            val lines = Source.fromFile(files.head).getLines().toList
            val predicted = lines(5).drop(20)
            
            if (predicted == "NONE"){
              pw.println(study.getName + "\t" + s.getName() + "\t" + predicted)
            } else {
              val lineages = predicted.split(", ").map(_.split("\\(").head)
              val covPerLin = lines.drop(8).grouped(4).map{linInfo =>
                (linInfo(0).replaceAll("\\s", "").split("->").last -> linInfo(1).replaceAll("\\s", "").split(":").last)
              }.toMap
              
              //println(predicted)
              val freqPerLin = predicted.split(", ").map{l => val arr = l.split("\\("); (arr(0) -> arr(1).dropRight(1))}.toMap
              //println(freqPerLin)
              
              val predicted2 = lineages.map(lin => lin + "(" + freqPerLin(lin) + "," + covPerLin(lin) + ")").mkString(", ")
              pw.println(study.getName + "\t" + s.getName() + "\t" + predicted2)
            }
            
          }
        }
      }

      pw.close
      
    }
  }
}