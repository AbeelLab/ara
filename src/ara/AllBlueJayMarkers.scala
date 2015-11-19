package ara

import java.io.File
import java.io.PrintWriter
import scala.io.Source

object AllBlueJayMarkers {
  
  case class Config(val dir: File = null, val output: String = null)
  
  def main (args: Array[String]) {
    val parser = new scopt.OptionParser[Config]("java -jar ara.jar all-bj-markers") {
      opt[File]('d', "dir") required () action { (x, c) => c.copy(dir = x) } text ("Directory with all BlueJay markers files.")
      opt[File]('o', "output") required () action { (x, c) => c.copy(output = x + ".txt") } text ("Output prefix.")
    }

    parser.parse(args, Config()) map { config =>
      
      val markerFiles = config.dir.listFiles().filter(_.getName.endsWith("_bluejaymarkers"))
      println(markerFiles.size + " marker files\n")
      
      val pw = new PrintWriter(config.output)      
      pw.println("# All markers for each hierarchical clusters.")
      pw.println("# Command: java -jar ara.jar all-bj-markers " + args.mkString(" "))
      
      var count = 0
      markerFiles.foreach{ file =>
        println("Processing " + file.getName + "...")
        val it = Source.fromFile(file).getLines.filterNot(_.startsWith("#")).grouped(2)        
        var fileCount = 0
        it.foreach{m => pw.println(m.mkString("\n")); fileCount = fileCount + 1}        
        println("\t" + fileCount + " markers")
        count = count + fileCount
      }
      println(count + " total markers")
      pw.close
    }
    
  }
}