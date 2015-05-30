//package macaw2
import java.io.File
import java.io.PrintWriter

/**
 * Script to write paths to VCFs in a text file. The output may be used for the script "vcf2snpAlignment.pl".
 */
object WriteVCFpaths {
  val usage = "usage: scala WriteVCFpaths.scala [directory] [output.txt]" 
  
  def main(args: Array[String]) {    
    def printToFile(f: File)(op: PrintWriter => Unit) {
      val writer = new PrintWriter(f)
      try { op(writer) } finally { writer.close() }
    }
    def listFiles(f: Any): List[File] = f match {
      case f: File if (f.isDirectory()) => f.listFiles().toList.flatMap(listFiles(_))
      case f: File if (f.isFile() && f.getName.endsWith(".vcf")) => List(f)
      case _ => Nil
    }
    
    args.length match {
      case 2 => {
        val list = listFiles(new File(args(0)))
        printToFile(new File(args(1))) { p => list.foreach(f => p.println(f)) }
      }
      case _ => println(usage)
    }
  }
}