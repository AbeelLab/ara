package ara

import java.io.File
import java.io.PrintWriter
import scala.io.Source

object BamInfo {
  val usage = "scala BamInfo.scala [directory] [out.txt]"

  def listFiles(f: Any): List[File] = f match {
    case f: File if (f.isDirectory()) => f.listFiles().toList.flatMap(listFiles(_))
    case f: File if (f.isFile() && f.getName.equals("reduced.vcf")) => List(f)
    case _ => Nil
  }

  def main(args: Array[String]) {
    
    if (args.length < 2) {
      println(usage)
    }
    else {
      val list = (0 to args.length -2).toList.flatMap(idx => listFiles(new File(args(idx)))).sortBy(file => file.getParentFile.getName)
      
      val pw = new PrintWriter(args(args.length - 1))
      list.foreach{vcf =>
        val cmd = Source.fromFile(vcf).getLines.toList(3).split(" ").drop(8).toList.mkString(" ")
       // val fc = cmd.count(_.equals("--frags"))
       // val jc = cmd.count(_.equals("--jumps"))
       // println(cmd.size)
        pw.println(cmd)
      }
      
      
      
      pw.close
    }
    
  }
}