package ara;
//package macaw2
import java.io.File
import java.io.PrintWriter

/**
 * Script to write paths to VCFs in a text file. The output may be used for the script "vcf2snpAlignment.pl".
 * Only VCFs named "reduced.vcf" will be filtered from the given directories.
 */
 
object WriteVCFpaths {
  val usage = "usage: scala WriteVCFpaths.scala [directory] [output.txt] \nMultiple directories may be given as input and should be separated with spaces." 
  
  def main(args: Array[String]) {    
    def printToFile(f: File)(op: PrintWriter => Unit) {
      val writer = new PrintWriter(f)
      try { op(writer) } finally { writer.close() }
    }
    def listFiles(f: Any): List[File] = f match {
      case f: File if (f.isDirectory()) => f.listFiles().toList.flatMap(listFiles(_))
      case f: File if (f.isFile() && f.getName.equals("reduced.vcf")) => List(f)
      case _ => Nil
    }
    
    args.length match {
      case n if (n > 1) => {
        val list = (0 to n-2).toList.flatMap(idx => listFiles(new File(args(idx)))).sortBy(file => file.getParentFile.getName)
        printToFile(new File(args(n-1))) { p => list.foreach(file => p.println(file))          
        }
      }
      case _ => println(usage)
    }
  }
}